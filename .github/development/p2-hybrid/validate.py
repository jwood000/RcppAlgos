#!/usr/bin/env python3
"""Build/run the temporary P2 checker, separately from benchmark timing."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess


HERE = Path(__file__).resolve().parent
SHARED = ("PiMeissel.cpp", "PiLehmer.cpp", "PhiParallel.cpp", "P3.cpp",
          "Phi.cpp", "Eratosthenes.cpp")
DEFAULT_VARIANTS = ("baseline", "2", "4", "8")
VARIANTS = DEFAULT_VARIANTS + ("16", "32", "64", "128", "selector", "selector-v2", "selector-v3")
SANITIZERS = {"asan-ubsan": "address,undefined", "tsan": "thread"}


def run(command, log, env=None):
    result = subprocess.run(command, text=True, capture_output=True,
                            timeout=300, env=env)
    log.write_text(result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError(f"command failed ({result.returncode}): {command}\n"
                           f"{result.stdout}{result.stderr}")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_dir", type=Path)
    parser.add_argument("--phase", choices=("build", "run", "all"), default="all")
    parser.add_argument("--variants", nargs="+", choices=VARIANTS,
                        default=DEFAULT_VARIANTS)
    parser.add_argument("--stream-window", action="store_true",
                        help="also check private stream windows through the maximum 128L endpoint")
    parser.add_argument("--selector-check", action="store_true",
                        help="also check the frozen selector policy and kernel forwarding")
    parser.add_argument("--selector-v2-check", action="store_true",
                        help="also check the frozen interval-aware v2 policy and kernel forwarding")
    parser.add_argument("--selector-v3-check", action="store_true",
                        help="also check the frozen v3 policy and kernel forwarding")
    args = parser.parse_args()
    if len(set(args.variants)) != len(args.variants):
        parser.error("--variants must not contain duplicates")
    root = args.artifact_dir.resolve()
    source = root / "source"
    output = root / "validation"
    output.mkdir(parents=True, exist_ok=True)
    if not (source / "src/Phi.cpp").is_file():
        parser.error("artifact directory must contain the immutable source snapshot")
    commands, results, stream_results, selector_results = [], [], [], []
    selector_v2_results = []
    selector_v3_results = []
    for name, sanitizer in SANITIZERS.items():
        directory = output / name
        directory.mkdir(exist_ok=True)
        flags = ["clang++", "-std=c++11", "-O1", "-g", "-pthread",
                 "-Wall", "-Wextra", "-Wpedantic", "-fno-omit-frame-pointer",
                 "-fno-sanitize-recover=all", f"-fsanitize={sanitizer}",
                 "-I" + str(source / "inst/include"), "-I" + str(HERE)]
        shared = []
        if args.phase != "run":
            units = [(source / "src" / filename, []) for filename in SHARED]
            units += [(source / "src/P2.cpp", ["-DP2=P2Reference"]),
                      (HERE / "check.cpp", [])]
            for path, extra in units:
                object_name = "P2Reference.o" if extra else path.stem + ".o"
                target = directory / object_name
                command = flags + extra + ["-c", str(path), "-o", str(target)]
                commands.append(command)
                run(command, directory / (object_name + ".log"))
                shared.append(target)
            for variant in args.variants:
                path = source / "src/P2.cpp" if variant == "baseline" else HERE / "P2Hybrid.cpp"
                extra = ([] if variant == "baseline" else ["-DP2_HYBRID_SELECTOR"]
                         if variant == "selector" else ["-DP2_HYBRID_SELECTOR_V2"]
                         if variant == "selector-v2" else ["-DP2_HYBRID_SELECTOR_V3"]
                         if variant == "selector-v3" else ["-DP2_HYBRID_MULTIPLIER=" + variant])
                target = directory / ("check-" + variant)
                command = flags + extra + [str(path)] + list(map(str, shared)) + ["-o", str(target)]
                commands.append(command)
                run(command, directory / ("build-" + variant + ".log"))
            print(f"Built {len(args.variants)} {name} checker binaries", flush=True)
            if args.stream_window:
                command = flags + ["-DP2_HYBRID_STREAM_CHECK",
                                   "-DP2_HYBRID_MULTIPLIER=128", str(HERE / "check.cpp")]
                command += [str(p) for p in shared if p.name != "check.o"]
                command += ["-o", str(directory / "check-stream-window")]
                commands.append(command)
                run(command, directory / "build-stream-window.log")
            if args.selector_check:
                command = flags + ["-DP2_HYBRID_SELECTOR", "-DP2_HYBRID_MEMORY",
                                   str(HERE / "selector-check.cpp"), str(HERE / "P2Hybrid.cpp")]
                command += [str(p) for p in shared if p.name != "check.o"]
                command += ["-o", str(directory / "check-selector-policy")]
                commands.append(command)
                run(command, directory / "build-selector-policy.log")
            if args.selector_v2_check:
                command = flags + ["-DP2_HYBRID_SELECTOR_V2", "-DP2_HYBRID_MEMORY",
                                   str(HERE / "selector-v2-check.cpp"), str(HERE / "P2Hybrid.cpp")]
                command += [str(p) for p in shared if p.name != "check.o"]
                command += ["-o", str(directory / "check-selector-v2-policy")]
                commands.append(command)
                run(command, directory / "build-selector-v2-policy.log")
            if args.selector_v3_check:
                command = flags + ["-DP2_HYBRID_SELECTOR_V3", "-DP2_HYBRID_MEMORY",
                                   str(HERE / "selector-v3-check.cpp"), str(HERE / "P2Hybrid.cpp")]
                command += [str(p) for p in shared if p.name != "check.o"]
                command += ["-o", str(directory / "check-selector-v3-policy")]
                commands.append(command)
                run(command, directory / "build-selector-v3-policy.log")

        if args.phase != "build":
            environment = os.environ.copy()
            environment.update(ASAN_OPTIONS="detect_leaks=0:halt_on_error=1",
                               UBSAN_OPTIONS="halt_on_error=1:print_stacktrace=1",
                               TSAN_OPTIONS="halt_on_error=1")
            for variant in args.variants:
                command = [str(directory / ("check-" + variant))]
                commands.append(command)
                result = run(command, directory / ("run-" + variant + ".log"), environment)
                matches = re.findall(r"^PASS checks=(\d+)$", result.stdout, re.MULTILINE)
                if len(matches) != 1 or result.stderr:
                    raise RuntimeError("unexpected checker output: " + result.stdout + result.stderr)
                results.append(dict(sanitizer=name, variant=variant, checks=int(matches[0])))
                print(f"{name} {variant}: {matches[0]} checks passed", flush=True)
            if args.stream_window:
                command = [str(directory / "check-stream-window")]
                commands.append(command)
                result = run(command, directory / "run-stream-window.log", environment)
                matches = re.findall(r"^PASS stream_checks=(\d+)$", result.stdout, re.MULTILINE)
                if len(matches) != 1 or result.stderr:
                    raise RuntimeError("unexpected stream checker output: " + result.stdout + result.stderr)
                stream_results.append(dict(sanitizer=name, checks=int(matches[0])))
                print(f"{name} stream windows: {matches[0]} checks passed", flush=True)
            if args.selector_check:
                command = [str(directory / "check-selector-policy")]
                commands.append(command)
                result = run(command, directory / "run-selector-policy.log", environment)
                matches = re.findall(r"^PASS selector_checks=(\d+) choice_cases=(\d+) kernel_cases=(\d+)$",
                                     result.stdout, re.MULTILINE)
                if len(matches) != 1 or result.stderr:
                    raise RuntimeError("unexpected selector checker output: " + result.stdout + result.stderr)
                count, choice_cases, kernel_cases = map(int, matches[0])
                selector_results.append(dict(sanitizer=name, checks=count,
                                             choice_cases=choice_cases, kernel_cases=kernel_cases))
                print(f"{name} selector policy: {count} checks passed", flush=True)
            if args.selector_v2_check:
                command = [str(directory / "check-selector-v2-policy")]
                commands.append(command)
                result = run(command, directory / "run-selector-v2-policy.log", environment)
                matches = re.findall(r"^PASS selector_v2_checks=(\d+) choice_cases=(\d+) kernel_cases=(\d+)$",
                                     result.stdout, re.MULTILINE)
                if len(matches) != 1 or result.stderr:
                    raise RuntimeError("unexpected v2 checker output: " + result.stdout + result.stderr)
                count, choice_cases, kernel_cases = map(int, matches[0])
                selector_v2_results.append(dict(sanitizer=name, checks=count,
                                                choice_cases=choice_cases, kernel_cases=kernel_cases))
                print(f"{name} selector v2 policy: {count} checks passed", flush=True)
            if args.selector_v3_check:
                command = [str(directory / "check-selector-v3-policy")]
                commands.append(command)
                result = run(command, directory / "run-selector-v3-policy.log", environment)
                matches = re.findall(r"^PASS selector_v3_checks=(\d+) choice_cases=(\d+) kernel_cases=(\d+)$",
                                     result.stdout, re.MULTILINE)
                if len(matches) != 1 or result.stderr:
                    raise RuntimeError("unexpected v3 checker output: " + result.stdout + result.stderr)
                count, choice_cases, kernel_cases = map(int, matches[0])
                selector_v3_results.append(dict(sanitizer=name, checks=count,
                                                choice_cases=choice_cases, kernel_cases=kernel_cases))
                print(f"{name} selector v3 policy: {count} checks passed", flush=True)

    (output / ("commands-" + args.phase + ".json")).write_text(json.dumps(commands, indent=2) + "\n")
    if results:
        if len({row["checks"] for row in results}) != 1:
            raise RuntimeError("checker paths differ between builds")
        (output / "results.json").write_text(json.dumps(results, indent=2) + "\n")
    if stream_results:
        (output / "stream-results.json").write_text(json.dumps(stream_results, indent=2) + "\n")
    if selector_results:
        (output / "selector-results.json").write_text(json.dumps(selector_results, indent=2) + "\n")
    if selector_v2_results:
        (output / "selector-v2-results.json").write_text(json.dumps(selector_v2_results, indent=2) + "\n")
    if selector_v3_results:
        (output / "selector-v3-results.json").write_text(json.dumps(selector_v3_results, indent=2) + "\n")
    paths = [HERE / "check.cpp", HERE / "P2Hybrid.cpp", HERE / "validate.py"]
    if "selector" in args.variants or args.selector_check:
        paths += [HERE / "P2Selector.h", HERE / "selector-check.cpp"]
    if "selector-v2" in args.variants or args.selector_v2_check:
        paths += [HERE / "P2SelectorV2.h", HERE / "selector-v2-check.cpp"]
    if "selector-v3" in args.variants or args.selector_v3_check:
        paths += [HERE / "P2SelectorV3.h", HERE / "selector-v3-check.cpp"]
    paths += [source / "src" / filename for filename in SHARED + ("P2.cpp",)]
    (output / "sources.json").write_text(json.dumps({str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                                                   for p in paths}, indent=2) + "\n")


if __name__ == "__main__":
    main()
