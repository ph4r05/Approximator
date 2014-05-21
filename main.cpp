/* 
 * File:   main.cpp
 * Author: ph4r05
 *
 * Created on May 16, 2014, 10:20 AM
 */

#include <cstdlib>

#include "base.h"
#include "Approximation.h"
#include "AESCipher.h"
#include "CombinatiorialGenerator.h"

#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

namespace po = boost::program_options;

using namespace std;
using namespace boost;

/**
 * Main entry point
 */
int main(int argc, char** argv) {
    AESCipher c;
    Approximation ap;
    
    po::options_description description("High order approximation");
    description.add_options()
            ("version,v",                                                               "Display the version number")
            ("help,h",                                                                  "Display this help message")
            ("corr,c",         po::value<ulong>()->default_value(0)->implicit_value(0), "Testing approximation correctness")
            ("acc",            po::value<ulong>()->default_value(0)->implicit_value(0), "Testing approximation accuracy samples")
            ("order,o",        po::value<uint>()->default_value(3)->implicit_value(3),  "Order limit for approximation")
//            ("extEnc,e",       po::value<bool>()->default_value(false)->implicit_value(false), "Use external encoding?")
//            ("out-file,o",     po::value<std::string>(),                                       "Output file to write encrypted data")
//            ("input-files",    po::value<std::vector<std::string>>(),                          "Input files")
//            ("create-table",   po::value<std::string>(),                                       "Create encryption/decryption tables");
    ;
            
    
    po::positional_options_description p;
    //p.add("input-files", -1);

    po::variables_map vm;
    po::store(po::command_line_parser(argc, argv).options(description).positional(p).run(), vm);
    po::notify(vm);

    if(vm.count("help")){
        std::cout << description;
        return 0;
    }

    if(vm.count("version")){
        std::cout << "High order approximation v0.1" << std::endl;
        return 0;
    }
    
    uint orderLimit = vm["order"].as<uint>();
    if (orderLimit > 5){
        cerr << "Approximation limit is too high, cannot continue" << endl;
        return 1;
    }
    
    // Prepare approximation object.
    ap.setCipher(&c);
    ap.setOrderLimit(orderLimit);
    ap.init();
    
    // Compute coefficients.
    ap.computeCoefficients();
        
    // Test approximation correctness
    ulong corrSamples = vm["corr"].as<ulong>();
    if (corrSamples>0){
        cout << "Testing approximation correctness: " << endl << " ";
        ap.selftestApproximation(corrSamples);
    }
    
    // Test approximation quality
    ulong accSamples = vm["acc"].as<ulong>();
    if (accSamples>0){
        // Here we test the accuracy of the high order approximation, random key,
        // random message. 
        cout << "Testing approximation quality: " << endl << " ";
        ap.testPolynomialApproximation(accSamples);
    }
    
    // Grober basis stuff.
    
    
    return 0;
}

