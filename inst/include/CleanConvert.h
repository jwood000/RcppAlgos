#ifndef CLEAN_CONVERT_H
#define CLEAN_CONVERT_H

#include <importExportMPZ.h>

namespace CleanConvert {

    template <typename stdType>
    void convertPrimitive(SEXP input, stdType &result,
                          std::string myErrorMsg) {
        
        switch(TYPEOF(input)) {
            case REALSXP: {
                result = Rcpp::as<stdType>(input);
                break;
            }
            case INTSXP: {
                result = Rcpp::as<stdType>(input);
                break;
            }
            case RAWSXP:
            case STRSXP: {
                mpz_t temp[1];
                mpz_init(temp[0]);
                createMPZArray(input, temp, 1);
                double dblTemp = mpz_get_d(temp[0]);
                result = dblTemp;
                mpz_clear(temp[0]);
                break;
            }
            default:
                Rcpp::stop(myErrorMsg);
        }
    }
    
    template <typename stdType>
    void convertVector(SEXP input, std::vector<stdType> &result,
                       std::string myErrorMsg) {
        
        int total = Rf_length(input);
        
        switch(TYPEOF(input)) {
            case REALSXP: {
                result = Rcpp::as<std::vector<stdType> >(input);
                break;
            }
            case INTSXP: {
                result = Rcpp::as<std::vector<stdType> >(input);
                break;
            }
            case RAWSXP: {
                const char* raw = (char*)RAW(input);
                total = ((int*)raw)[0];
                // do not put a break here. Fall to
                // the next case for complete conversion
            }
            case STRSXP: {
                mpz_t *temp;
                temp = (mpz_t *) malloc(sizeof(mpz_t) * total);
                
                for (int i = 0; i < total; ++i)
                    mpz_init(temp[i]);
                
                createMPZArray(input, temp, total);
                std::vector<double> dblTemp(total);
                
                for (int i = 0; i < total; ++i) {
                    dblTemp[i] = mpz_get_d(temp[i]);
                    result.push_back((stdType) dblTemp[i]);
                }
                
                for (int i = 0; i < total; ++i)
                    mpz_clear(temp[i]);
                
                free(temp);
                break;
            }
            default:
                Rcpp::stop(myErrorMsg);
        }
    }
}

#endif
