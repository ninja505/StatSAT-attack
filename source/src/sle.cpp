#include "sle.h"
#include "sld.h"
#include "encoder.h"
#include "mutability.h"
#include "dac12enc.h"
#include "randomins.h"
#include "optenc.h"
#include "toc13enc.h"
#include <math.h>

#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <stdio.h>
#include <signal.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <vector>
#include <regex>
#include <boost/algorithm/string/replace.hpp>

std::string output_file;
int target_keys = -1;

// typical usages: 
//    -M <enc> to compute mutability graph and dump it to a file.
//    -d -g <graph> -C  <clique> -o <enc> <bench>
//    -r 1 -k <keys> to do random insertion.
//    -r 1 -f <fraction> to do random insertion.
//    -I <ILP> encoder -f <fraction> -o <encoded-file> <bench>
//    -t -f <fraction> <bench> to dump the fault impact files.
//    -t -f <fraction> -s <bench> to dump the fault impact files and for mux encoding
//    -t -f <fraction> -T <faultimpact> -o <output> <bench> : to encode using fault impact.
//    -t -f <fraction> -s -T <faultimpact> -o <output> <bench> : to encode using fault impact with muxes.
//    -i -f <fraction> -i -o <output> <bench> : to encode using the IOLTS technique.
int sle_main(int argc, char* argv[])
{
    int cpu_limit = -1;
    int64_t data_limit = -1;
    bool extended = true;
    std::string mutability;
    std::string graph;
    std::string clique;
    int dac12_enc = 0;
    int random_ins = 0;
    int ilp_encoder = 0;
    int toc13_enc = 0;
    std::string fault_impact_file;
    double key_fraction = 0.0;
    int mux_enc = 0;
    int iolts14_enc = 0;
    bool error_prop = false;

    int c;
    while ((c = getopt (argc, argv, "ihpc:r:m:o:k:eM:dg:C:f:ItT:s")) != -1) {
        switch (c) {
            case 'h':
                return sle_usage(argv[0]);
                break;
            case 'i':
                iolts14_enc = 1;
                break;
            case 't':
                toc13_enc = 1;
                break;
            case 's':
                mux_enc = 1;
                break;
            case 'T':
                fault_impact_file = optarg;
                break;
            case 'c':
                cpu_limit = atoi(optarg);
                break;
            case 'm':
                data_limit = ((int64_t) atoi(optarg)) * ((int64_t)1048576);
                break;
            case 'o':
                output_file = optarg;
                break;
            case 'k':
                target_keys = atoi(optarg);
                break;
            case 'f':
                key_fraction = atof(optarg);
                break;
            case 'M':
                mutability = optarg;
                break;
            case 'g':
                graph = optarg;
                break;
            case 'C':
                clique = optarg;
                break;
            case 'e':
                extended = !extended;
                break;
            case 'd':
                dac12_enc = 1;
                break;
            case 'I':
                ilp_encoder = 1;
                break;
            case 'r':
                random_ins = atoi(optarg);
                break;
            case 'p':
                error_prop = true;
            default:
                break;
        }
    }


    // check if we got a test article.
    if (error_prop) {
        if (optind == argc || optind != argc-3) 
            return sle_usage(argv[0]);
    }
    else if(optind == argc || optind != argc-1) {
        return sle_usage(argv[0]);
    }

    if(cpu_limit != -1) {
        setup_cpux_handler();

        struct rlimit rl;
        rl.rlim_cur = cpu_limit;
        setrlimit(RLIMIT_CPU, &rl);
    }

    if(data_limit != -1) {
        struct rlimit rl;
        rl.rlim_cur = rl.rlim_max = data_limit;
        if(setrlimit(RLIMIT_AS, &rl) != 0) {
            std::cerr << "trouble setting data limit." << std::endl;
        }
    }

    yyin = fopen(argv[optind], "rt");
    if(yyin == NULL) {
        perror(argv[optind]);
        return 1;
    }

    if(yyparse() == 0) {
        using namespace ast_n;
    
        if(iolts14_enc) {
            if(key_fraction == 0.0) {
                std::cerr << "Error: must specify fraction to insert. " << std::endl;
                exit(1);
            }
            ckt_n::toc13enc_t tenc(*statements, key_fraction);
            tenc.encodeIOLTS14();
            if(output_file.size() == 0) {
                tenc.write(std::cout);
            } else {
                std::ofstream fout(output_file.c_str());
                tenc.write(fout);
            }
        } else if(toc13_enc) {
            if(key_fraction == 0.0) {
                std::cerr << "Error: must specify fraction to insert. " << std::endl;
                exit(1);
            }
            ckt_n::toc13enc_t tenc(*statements, key_fraction);
            if(fault_impact_file.size() == 0) {
                tenc.evaluateFaultImpact(5000);
            } else {
                tenc.readFaultImpact(fault_impact_file);
            }
            if(mux_enc) {
                tenc.encodeMuxes();
            } else {
                tenc.encodeXORs();
            }
            if(output_file.size() == 0) {
                tenc.write(std::cout);
            } else {
                std::ofstream fout(output_file.c_str());
                tenc.write(fout);
            }

        } else if(ilp_encoder) {
            if(key_fraction == 0.0) {
                std::cerr << "Error: must specify fraction to insert. " << std::endl;
                exit(1);
            }
            ckt_n::ilp_encoder_t ioenc(*statements, key_fraction);
            ioenc.encode3();
            if(output_file.size() == 0) {
                ioenc.write(std::cout);
            } else {
                std::ofstream fout(output_file.c_str());
                ioenc.write(fout);
            }
        } else if(dac12_enc) {
            if(graph.size() == 0) {
                std::cerr << "Error: must specify graph filename." << std::endl;
                exit(1);
            }
            if(clique.size() == 0) {
                std::cerr << "Error: must specify clique filename." << std::endl;
                exit(1);
            }
            if(output_file.size() == 0) {
                std::cerr << "Error: must specify output file." << std::endl;
                exit(1);
            }
            if(key_fraction == 0.0) {
                std::cerr << "Error: must specify fraction to insert. " << std::endl;
                exit(1);
            }
            ckt_n::dac12enc_t denc(*statements, graph, clique, key_fraction);
            std::ofstream fout(output_file.c_str());
            denc.write(fout);
        } else if(mutability.size()) {
            ckt_n::ckt_t ckt(*statements);

            std::ofstream fout(mutability.c_str());
            ckt_n::mutability_analysis_t mut(ckt, fout);
            mut.analyze();
            fout.close();

            std::cout << "finished" << std::endl;
        } else if(random_ins == 1) {
            if(target_keys == -1 && key_fraction == 0.0) {
                std::cerr << "Error, must specify target number of keys with -k <keys> or -f <fraction> flag." << std::endl;
            ckt_n::ckt_t ckt(*statements);

            if(key_fraction != 0.0) {
                target_keys = (int) (ckt.num_gates() * key_fraction + 0.5);
            }
            random_insert(ckt, target_keys);
            std::cout << ckt << std::endl;
          }
        } else if(error_prop == true) {
            ckt_n::ckt_t ckt(*statements);
            // probability and mapping files
            std::fstream fprob;
            std::ifstream ffmap; 
            std::ofstream f_tmp;
            // maps of flip flop d and q pins
            std::map<std::string, std::string> ff_dmap;
            std::map<std::string, std::string> ff_qmap;
            // map for output name to output nodelist index
            std::map<std::string, int> output_idx_map;
            // signal and error probability vectors
            std::vector<double> v_sig_prob(ckt.num_inputs(), 0.5);
            std::vector<double> v_err_prob(ckt.num_inputs(), 0);
            std::pair<std::vector<double>, std::vector<double>> output_probs;
            std::vector<double> v_out_sig;
            std::vector<double> v_out_err;
            std::string line;

            // parse ffmap to create flip flop i/o pin mapping
            ffmap.open(argv[optind+2]);
            if (ffmap.is_open())
            {
                while (std::getline(ffmap, line))
                {
                    std::istringstream csvline(line); 
                    std::vector<std::string> tokens;
                    std::string str; 
                    while (std::getline(csvline, str, ','))
                    {
                        tokens.push_back(str);
                    }
                    ff_dmap[tokens[0]] = tokens[1];
                    ff_qmap[tokens[0]] = tokens[2];
                }
                ffmap.close();
            }

            // parse probability file for inputs
            fprob.open(argv[optind+1]);
            if (fprob.is_open())
            {
                while (std::getline(fprob, line))
                {
                    double sig_prob, err_prob;
                    std::string name;

                    std::istringstream csvline(line);
                    std::vector<std::string> tokens;
                    std::string str;
                    std::regex re("[^/]*$");
                    std::regex node_re_begin("[^f].*");
                    std::regex node_re_end("_[0-9]+_$");
                    std::smatch match;
                    while (std::getline(csvline, str, ','))
                    {
                        tokens.push_back(str);
                    }
                    std::regex_search(tokens[0], match, re);

                    name = match[0];
                    sig_prob = std::stod(tokens[1]);
                    err_prob = std::stod(tokens[2]);
                    
                    // different cases for inputs and flip flops
                    // ignore outputs
                    if (tokens[0].find("ports_in") != std::string::npos)
                    {
                        boost::replace_all(name, "[", "0");
                        boost::replace_all(name, "]", "0");
                        for (int i = 0; i < ckt.num_inputs(); i++)
                        {
                            std::string node_name = ckt.inputs[i]->name;
                            
                            // node regex match: take away "f" at beginning and _[0-9]+_ at end
                            std::regex_search(node_name, match, node_re_begin);
                            node_name = match[0];
                            std::regex_search(node_name, match, node_re_end);
                            str = match[0];
                            boost::replace_all(node_name, str, "");
                            
                            if (node_name == name)
                            {
                                //std::cout << name + " matches " + node_name << std::endl;
                                v_sig_prob[i] = sig_prob;
                                v_err_prob[i] = err_prob;
                            }
                        }
                    }
                    else if (tokens[0].find("instances_seq") != std::string::npos)
                    {
                        // q pins are circuit inputs
                        std::string in_name = ff_qmap[name]; 
                        boost::replace_all(in_name, "[", "0");
                        boost::replace_all(in_name, "]", "0");
                        for (int i = 0; i < ckt.num_inputs(); i++)
                        {
                            std::string node_name = ckt.inputs[i]->name;
                            
                            // find a better way for dealing with the o's
                            boost::replace_all(node_name, "Interfaceo", "Interface");
                            
                            // node regex match: take away "f" at beginning and _[0-9]+_ at end
                            std::regex_search(node_name, match, node_re_begin);
                            node_name = match[0];
                            std::regex_search(node_name, match, node_re_end);
                            str = match[0];
                            boost::replace_all(node_name, str, "");

                            if (node_name == in_name)
                            {
                                //std::cout << in_name + " matches " + node_name << std::endl;
                                v_sig_prob[i] = sig_prob;
                                v_err_prob[i] = err_prob;
                            }
                        }

                        // d pins are circuit outputs
                        std::string out_name = ff_dmap[name];
                        boost::replace_all(out_name, "[", "0");
                        boost::replace_all(out_name, "]", "0");
                        for (int i = 0; i < ckt.num_outputs(); i++)
                        {
                            std::string node_name = ckt.outputs[i]->name;

                            // node regex match: take away "f" at beginning and _[0-9]+_ at end
                            std::regex_search(node_name, match, node_re_begin);
                            node_name = match[0];
                            std::regex_search(node_name, match, node_re_end);
                            str = match[0];
                            boost::replace_all(node_name, str, "");

                            if (node_name == out_name)
                            {
                                //std::cout << out_name + " matches " + node_name << std::endl;
                                // note that name, not out_name, is used here;
                                // this is so only one map needs to be indexed when
                                // writing to output
                                output_idx_map[name] = i;
                            }
                        }
                    }
                    else // ports_out
                    {
                        boost::replace_all(name, "[", "0");
                        boost::replace_all(name, "]", "0");

                        for (int i = 0; i < ckt.num_outputs(); i++)
                        {
                            std::string node_name = ckt.outputs[i]->name;
                            // node regex match: take away "f" at beginning and _[0-9]+_ at end
                            std::regex_search(node_name, match, node_re_begin);
                            node_name = match[0];
                            std::regex_search(node_name, match, node_re_end);
                            str = match[0];
                            boost::replace_all(node_name, str, "");

                            if (node_name == name)
                            {
                                //std::cout << name + " matches " + node_name << std::endl;
                                output_idx_map[name] = i;
                            }
                        }
                    }
                }
                fprob.close();
            }
            
            // call evaluate_probs_as to get output signal and error probabilities 
            output_probs = ckt.evaluate_probs_as(v_sig_prob, v_err_prob);
            v_out_sig = output_probs.first;
            v_out_err = output_probs.second;

            // print output to file
            fprob.open(argv[optind+1]);
            f_tmp.open("error_prop.tmp");
            if (fprob.is_open() && f_tmp.is_open())
            {
                while (std::getline(fprob, line))
                {
                    std::string name;
                    std::istringstream csvline(line);
                    std::vector<std::string> tokens;
                    std::string str;
                    std::regex re("[^/]*$");
                    std::smatch match;
                    int idx;

                    while (std::getline(csvline, str, ','))
                    {
                        tokens.push_back(str);
                    }

                    if (tokens[0].find("ports_out") != std::string::npos || 
                        tokens[0].find("instances_seq") != std::string::npos)
                    {
                        std::regex_search(tokens[0], match, re);
                        name = match[0];
                        idx = output_idx_map[name]; 

                        // write line with updated probs 
                        f_tmp << tokens[0] + ',' + std::to_string(v_out_sig[idx]) + ',' + \
                            std::to_string(v_out_err[idx]) << std::endl;
                    }
                    else
                    {
                        // write same line
                        f_tmp << tokens[0] + ',' + tokens[1] + ',' + tokens[2] << std::endl; 
                    }

                }
                fprob.close();
                f_tmp.close();
            }
            // remove orig prob file
            std::remove(argv[optind+1]);
            // rename error_prop.tmp to orig prob file's name
            std::rename("error_prop.tmp", argv[optind+1]);
        } else {
            std::cout << "Error! Must select an encoding scheme. " << std::endl;
            exit(1);
        }
    }


    return 0;
}

int sle_usage(const char* progname)
{
    std::cout << "Usage: " << progname << " [options] <bench-file> [<signal-error-prob-file> <flip-flop-mapping-file>]" 
              << std::endl;
    std::cout << "Options may be one of the following." << std::endl;
    std::cout << "    -h            : this message." << std::endl;
    std::cout << "    -e            : toggle extended encoder (default=true)." << std::endl;
    std::cout << "    -o <filename> : output file." << std::endl;
    std::cout << "    -k <keys>     : number of keys to introduces (default=10% of num_gates)." << std::endl;
    std::cout << "    -c <value>    : CPU time limit (s)." << std::endl;
    std::cout << "    -m <value>    : mem usage limit (MB)." << std::endl;
    std::cout << "    -m <value>    : mem usage limit (MB)." << std::endl;
    std::cout << "    -p            : propagate error through circuit (requires extra files)." << std::endl;
    return 0;
}


// end of the file.
// 
//
//
