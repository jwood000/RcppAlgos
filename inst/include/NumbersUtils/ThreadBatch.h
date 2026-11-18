#pragma once

#include <cstddef>
#include <exception>
#include <stdexcept>
#include <thread>
#include <vector>

namespace PrimeCounting {

    // One batch owns its workers and their exception slots. Only the calling
    // thread may launch or finish the batch. Each copyable callable supplies
    // its own work; the batch knows nothing about Phi, tables, or result types.
    // Any objects captured by reference must outlive the batch, including
    // when a later thread fails to start. Declare result storage before it.
    class ThreadBatch {
        std::vector<std::exception_ptr> errors;
        std::vector<std::thread> threads;

        void join() {
            for (auto& thread : threads) {
                if (thread.joinable()) thread.join();
            }
        }

    public:
        explicit ThreadBatch(std::size_t capacity) : errors(capacity) {
            threads.reserve(capacity);
        }

        ThreadBatch(const ThreadBatch&) = delete;
        ThreadBatch& operator=(const ThreadBatch&) = delete;

        ~ThreadBatch() {
            // Also runs if a launch fails after earlier workers have started.
            join();
        }

        template <typename Function>
        void launch(Function work) {
            const std::size_t slot = threads.size();

            if (slot == errors.size()) {
                throw std::length_error("ThreadBatch capacity exceeded");
            }

            threads.emplace_back([this, work, slot]() mutable {
                try {
                    work();
                } catch (...) {
                    // Each worker writes only its fixed slot. Callables must
                    // not use the R API; propagate C++ exceptions instead.
                    errors[slot] = std::current_exception();
                }
            });
        }

        void finish() {
            join();
            // Joining makes worker writes visible. Propagate in launch order
            // on the calling thread, before the caller uses partial results.
            for (const auto& error : errors) {
                if (error) std::rethrow_exception(error);
            }
        }
    };
}
