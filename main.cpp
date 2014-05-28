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
    po::options_description description("High order approximation");
    description.add_options()
            ("version,v",                                                               "Display the version number.")
            ("help,h",                                                                  "Display this help message.")
            ("corr,c",         po::value<ulong>()->default_value(0)->implicit_value(0), "Testing approximation correctness.")
            ("acc",            po::value<ulong>()->default_value(0)->implicit_value(0), "Testing approximation accuracy samples.")
            ("order,o",        po::value<uint>()->default_value(3)->implicit_value(3),  "Order limit for approximation.")
            ("key,k",          po::value<bool>()->default_value(false)->implicit_value(false), "Try to solve key equations with GB.")
            ("itest",          po::value<bool>()->default_value(false)->implicit_value(false), "Internal implementation correctness tests.")
            ("key-to-zero,z",  po::value<uint>()->default_value(0)->implicit_value(0),  "Number of key-bits set to zero in GB evaluation.")
            ("polymap,p",      po::value<std::vector<std::string>>(),                   "List of polynomials to take into solution.")
            ("samples,s",      po::value<uint>()->default_value(1)->implicit_value(1),  "Number of samples to try.")
            ("dump-input",     po::value<bool>()->default_value(false)->implicit_value(false), "Whether to dump input base to GB.")
            ("self-test",      po::value<bool>()->default_value(false)->implicit_value(false), "Self-testing during normal computation.")
            ("basis-reduction",po::value<uint>()->default_value(0)->implicit_value(0), "Reduce over-determined systems to exactly determined.")
            ("rounds,r",       po::value<int>()->default_value(-1)->implicit_value(-1), "Number of rounds of the cipher.")
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
    
    bool selftest = vm["self-test"].as<bool>();
    
    // Prepare approximation object.
    AESCipher c;
    if (vm.count("rounds")){
        c.setNumRounds(vm["rounds"].as<int>());
    }
    
    // Approximation object - main one.
    Approximation ap(orderLimit);
    
    // Set cipher
    ap.setCipher(&c);
    
    // Number of key bits set to 0.
    ap.setKeybitsToZero(vm["key-to-zero"].as<uint>());
    
    ap.init();
    
    // Internal implementation test.
    bool itest = vm["itest"].as<bool>();
    if (itest){
        cout << "Testing internal implementation" << endl;
        ap.selftestIndexing();
        return 0;
    }
    
    // Compute coefficients.
    ap.computeCoefficients(ap.getCoefficients());
        
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
    
    // GB stuff.
    bool keyEq = vm["key"].as<bool>();
    if (keyEq){
        if(vm.count("polymap")){
            std::vector<std::string> polynomials2take = vm["polymap"].as<std::vector<std::string>>();
            cout << "Polymap set, size="<<polynomials2take.size()<<endl;
            ap.setPoly2Take(polynomials2take);
        }
        
        cout << "Solving key equations with GB" << endl;
        ap.initFGb(ap.getNumVariables());
        ap.solveKeyGrobner(vm["samples"].as<uint>(), vm["dump-input"].as<bool>(), selftest, vm["basis-reduction"].as<uint>()); 
        ap.deinitFGb();
    }
    
    return 0;
}

