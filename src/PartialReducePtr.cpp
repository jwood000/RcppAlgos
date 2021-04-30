#include <string>

template <typename T>
using partialReducePtr = void (*const)(int m, T &partial, T w);

template <typename T>
void PartialReduceProd(int m, T &partial, T w) {
    partial /= w;
}

template <typename T>
void PartialReduceSum(int m, T &partial, T w) {
    partial -= w;
}

template <typename T>
void PartialReduceMean(int m, T& partial, T w) {
    partial = (partial * static_cast<double>(m) - w) / static_cast<double>(m - 1);
}

template <typename T>
partialReducePtr<T> GetPartialReducePtr(const std::string &myFun) {
    
    if (myFun == "prod") {
        return(partialReducePtr<T>(PartialReduceProd));
    } else if (myFun == "sum") {
        return(partialReducePtr<T>(PartialReduceSum));
    } else {
        return(partialReducePtr<T>(PartialReduceMean));
    }
}

template partialReducePtr<int> GetPartialReducePtr(const std::string&);
template partialReducePtr<double> GetPartialReducePtr(const std::string&);
