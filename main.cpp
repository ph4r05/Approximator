/* 
 * File:   main.cpp
 * Author: ph4r05
 *
 * Created on May 16, 2014, 10:20 AM
 */

#include <cstdlib>
#include <fstream>
#include <iomanip>

#include "base.h"
#include "Approximation.h"
#include "AESCipher.h"
#include "CombinatiorialGenerator.h"
#include "Keccak2.h"
#include "KeccakOptAsm.h"
#include "Keccak1600.h"

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
            ("threads,t",      po::value<uint>()->default_value(1)->implicit_value(1),  "Number of threads to use for computation.")
            ("cube",           po::value<uint>()->default_value(0)->implicit_value(0),  "Starts cube attack.")
            ("alg",            po::value<uint>()->default_value(0)->implicit_value(0),  "Algorithm to analyze. 0=AES, 1=Keccak, 2=Keccak 1600 (no key), 3=Keccak assembler optimized")
            ("relations",      po::value<uint>()->default_value(128)->implicit_value(128),  "Number of relations finding rounds.")
            ("subcube",        po::value<uint>()->default_value(0)->implicit_value(0),  "Subcubes to compute in parallel.")
            ("wkey",           po::value<uint>()->default_value(1)->implicit_value(1),  "Weight of the key cube.")
            ("verbose",        po::value<uint>()->default_value(1)->implicit_value(1),  "Verbosivity level.")
            ("saverel",        po::value<bool>()->default_value(true)->implicit_value(true),  "Save relations to a file.")
            ("dumprel",        po::value<bool>()->default_value(false)->implicit_value(false),  "Dump all relations found (for statistical testing).")
            ("dump-approx",    po::value<std::string>(),                                "File where to dump approximation of the function.")
            ("dump-format",    po::value<uint>()->default_value(1)->implicit_value(1),  "Output dump format. 1=ASCII, 2=Binary, newline, 3=Binary") 
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
    
    // Algorithm to analyze.
    ICipher * c;
    
    uint algId = vm["alg"].as<uint>();
    switch(algId){
        case 0: 
            cout << "Algorithm=AES" << endl;
            c = new AESCipher();
            break;
        case 1:
            cout << "Algorithm=Keccak" << endl;
            c = new Keccak2();
            break;
        case 2:
            cout << "Algorithm=Keccak1600 (no key)" << endl;
            c = new Keccak1600();
            break;
        case 3:
            cout << "Algorithm=Keccak Assembler optimized" << endl;
            c = new KeccakOptAsm();
            break;
        default: 
            cerr << "Unknown algorithm id="<<algId<<endl;
            return 8;
    }
    
    // Prepare approximation object.
    if (vm.count("rounds")){
        int res = c->setNumRounds(vm["rounds"].as<int>());
        cout << "Setting number of rounds=" << (c->getNumRounds()) << endl;
        assert(res == 1);
    }
    
    // Seed random number generator
    unsigned seed = (unsigned)time(0);
    srand(seed); 
    cout << "Generator seeded with seed=" << seed << endl;
    
    // Approximation object - main one.
    Approximation ap(orderLimit);
    // verbosity level
    ap.setVerboseLvl(vm["verbose"].as<uint>());
    
    // Set cipher
    ap.setCipher(c);
    
    // Number of key bits set to 0.
    ap.setKeybitsToZero(vm["key-to-zero"].as<uint>());
    // Thread count to use.
    ap.setThreadCount(vm["threads"].as<uint>());
    // Internal approximator initialization.
    ap.init();
    
    // Cube attack switch.
    uint cube = vm["cube"].as<uint>();
    
    // Internal implementation test.
    bool itest = vm["itest"].as<bool>();
    if (itest){
        cout << "Testing internal implementation" << endl;
        ap.selftestIndexing();
        
        delete c;
        return 0;
    }
    
    // Compute coefficients.
    // If in cube attack mode, approximation is not needed.
    if (cube==0){
        cout << "Computing coefficients" << endl;
        ap.computeCoefficients(ap.getCoefficients());
        
        // If coefficient approximation is desired.
        if (vm.count("dump-approx")){
            std::string approxFile = vm["dump-approx"].as<std::string>();
            ofstream relOf(approxFile.c_str(), ios_base::app);
            ap.dumpOutputFunctions(relOf, ap.getCoefficients(), orderLimit, 
                    c->getKeyBlockSize()*8 + c->getInputBlockSize()*8,
                    c->getOutputBlockSize()*8,
                    false,
                    vm["dump-format"].as<uint>());
        }
    }
        
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
    
    // Cube stuff.
    if (cube>0){
        cout << "Cube attack" << endl;
        ap.initFGb(ap.getNumVariables());
        ap.cubeAttack(
                cube, 
                vm["wkey"].as<uint>(), 
                vm["relations"].as<uint>(), 
                vm["subcube"].as<uint>(), 
                vm["saverel"].as<bool>(),
                vm["dumprel"].as<bool>()
        );
        ap.deinitFGb();
    }
    
    delete c;
    return 0;
}

