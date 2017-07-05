/*
 * >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 * ART -- Artificial Read Transcription, ART_454
 * Authors: Weichun Huang 2008-2016
 * License: GPL v3 
 * ############################################################################
 * #    This program is free software: you can redistribute it and/or modify  #
 * #    it under the terms of the GNU General Public License as published by  #
 * #    the Free Software Foundation, either version 3 of the License, or     #
 * #    (at your option) any later version.                                   #
 * #                                                                          #
 * #    This program is distributed in the hope that it will be useful,       #
 * #    but WITHOUT ANY WARRANTY; without even the implied warranty of        #
 * #    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         #
 * #    GNU General Public License for more details.                          #
 * #                                                                          #
 * #    You should have received a copy of the GNU General Public License     #
 * #    along with this program.  If not, see <http://www.gnu.org/licenses/>. #
 * ############################################################################
 * <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*/

#include <iostream>
#include <sstream>
#include <string>
#include <time.h>
#include <algorithm>
#include <iomanip>
#include <ctime>
#include "art.h"
#include "samRead.h"
#include "readSeqFile.h"

using namespace std;

#define PRGNAME "art_454"
#define MINIMUM_FLAGMENT_SIZE 100

int main(int argc, char* argv[]){
    cout <<"==================================================================="<<endl;
    cout <<"                      ART_454 (Version 2.6.0)                      "<<endl; 
    cout <<"                   Simulation of 454 Pyrosequencing                "<<endl; 
    cout <<"     Copyright (c) 2008-2015, Weichun Huang. All Rights Reserved.  "<<endl; 
    cout <<"==================================================================="<<endl<<endl;

    bool amplicon=false;
    bool amplicon_pair=false;
    bool debug=false;
    bool fixed_seed = false;
    unsigned int rand_seed = 0;
    bool sam_out=false;
    bool aln_out=false;
    bool titanium=false;
    bool user_cycle=false;
    bool use_cigarM = false;
    string profile_name = "";
    int num_flow_cycles = 100;
    int i=1;
    for(;i<argc;++i){
	    char* pch = argv[i];
	    if( *pch != '-' || *(pch+1) == '\0') break;
	    while(*++pch){
		    switch(*pch){
			    case 'p':
			    case 'P':
				    if(i<argc) profile_name=argv[++i]; 
				    break;
			    case 'd':
				    debug=true;
				    break;
			    case 'a':
				    aln_out=true;
				    break;
			    case 'A':
				    amplicon=true;
				    break;
			    case 'B':
				    amplicon=true;
				    amplicon_pair=true;
				    break;
			    case 'M':
				    use_cigarM=true;
				    break;
			    case 's':
				    sam_out=true;
				    break;
			    case 't':
			    case 'T':
				    titanium=true;
				    break;
			    case 'c':
			    case 'C':
				    if(i<argc) num_flow_cycles = atoi(argv[++i]);
				    user_cycle=true;
				    break;
			    case 'r':
				    if(i<argc){
					    rand_seed =atoi(argv[++i]);
					    fixed_seed = true;
				    }
				    else {
					    cerr<<"Error: no random seed provided"<<endl;
					    exit(1);
				    }
				    break;
			    default:
				   cerr<<"Error: unreconized option \""<<*pch<< "\""<<endl;
				   break;
		    }
	    }
    }

    if(!user_cycle && titanium){
	   num_flow_cycles=200;
    }

    int k=argc-i+1;
    if(k !=4 &&  k !=6){
        cout<<"===== USAGE ===== "<<endl<<endl;
        cout<<"SINGLE-END SIMULATION"<<endl;
        cout<<"\t"<< PRGNAME <<" [-s] [-a ] [-t] [-r rand_seed] [ -p read_profile ] [ -c num_flow_cycles ] <INPUT_SEQ_FILE> <OUTPUT_FILE_PREFIX> <FOLD_COVERAGE>"<<endl<<endl;
        cout<<"PAIRED-END SIMULATION"<<endl;
        cout<<"\t"<< PRGNAME <<" [-s] [-a ] [-t] [-r rand_seed] [ -p read_profile ] [ -c num_flow_cycles ] <INPUT_SEQ_FILE> <OUTPUT_FILE_PREFIX> <FOLD_COVERAGE> <MEAN_FRAG_LEN> <STD_DEV>"<<endl<<endl;
        cout<<"AMPLICON SEQUENCING SIMULATION"<<endl;
        cout<<"\t"<< PRGNAME <<" [-s] [-a ] [-t] [-r rand_seed] [ -p read_profile ] [ -c num_flow_cycles ] <-A|-B> <INPUT_SEQ_FILE> <OUTPUT_FILE_PREFIX> <#_READS/#_READ_PAIRS_PER_AMPLICON>"<<endl<<endl;
        cout<<"===== PARAMETERS ====="<< endl<<endl;
        cout<<"INPUT_SEQ_FILE           -  the filename of DNA/RNA reference sequences in FASTA format"<<endl;
        cout<<"OUTPUT_FILE_PREFIX       -  the prefix or directory of output read data file (*.fq) and read alignment file (*.aln)"<<endl;
	cout<<"FOLD_COVERAGE            -  the fold of read coverage over the reference sequences"<<endl; 
	cout<<"MEAN_FRAG_LEN            -  the average DNA fragment size for paired-end read simulation"<<endl; 
	cout<<"STD_DEV                  -  the standard deviation of the DNA fragment size for paired-end read simulation"<<endl; 
        cout<<"#READS_PER_AMPLICON      -  number of reads per amplicon (for 5'end amplicon sequencing)"<<endl;
	cout<<"#READ_PAIRS_PER_AMPLICON -  number of read pairs per amplicon (for two-end amplicon sequencing)"<<endl<<endl; 
        cout<<"===== OPTIONAL PARAMETERS ====="<< endl<<endl;
        cout<<"-A indicate to perform single-end amplicon sequencing simulation"<<endl;
        cout<<"-B indicate to perform paired-end amplicon sequencing simulation"<<endl;
        cout<<"-M indicate to use CIGAR 'M' instead of '=/X' for alignment match/mismatch"<<endl;
        cout<<"-a indicate to output the ALN alignment file"<<endl;
        cout<<"-s indicate to output the SAM alignment file"<<endl;
        cout<<"-d print out warning messages for debugging"<<endl;
        cout<<"-t indicate to simulate reads from the built-in GS FLX Titanium profile [default: GS FLX profile]"<<endl;
        cout<<"-r specify a fixed random seed for the simulation (to generate two identical datasets from two different runs)"<<endl;
        cout<<"-c specify the number of flow cycles by the sequencer [ default: 100 for GS-FLX, and 200 for GS-FLX Titanium ] "<<endl;
        cout<<"-p specify user's own read profile for simulation"<<endl;
       	cout<<"   NOTE: the name of a read profile is the directory containing read profile data files."<<endl;
        cout<<"         please read the REAME file about the format of 454 read profile data files and."<<endl;
        cout<<"         and the default filenames of these data files." <<endl<<endl; 
        cout<<"===== EXAMPLES ====="<< endl<<endl;
        cout<<" 1) singl-end simulation with 20X coverage"<<endl;
       	cout<<"\t"<<PRGNAME <<" -s seq_reference.fa ./outdir/single_dat 20"<<endl;
        cout<<" 2) paired-end simulation with the mean fragment size 1500 and STD 20 using GS FLX Titanium platform"<<endl;
       	cout<<"\t"<<PRGNAME <<" -s -t seq_reference.fa ./outdir/paired_dat 10 1500 20"<<endl;
        cout<<" 3) paired-end simulation with a fixed random seed"<<endl;
       	cout<<"\t"<<PRGNAME <<" -s -r 777 seq_reference.fa ./outdir/paired_fxSeed 10 2500 50"<<endl;
        cout<<" 4) single-end amplicon sequencing with 10 reads per amplicon"<<endl;
       	cout<<"\t"<<PRGNAME <<" -A -s amplicon_ref.fa ./outdir/amp_single 10"<<endl;
        cout<<" 5) paired-end amplicon sequencing with 10 read pairs per amplicon"<<endl;
       	cout<<"\t"<<PRGNAME <<" -B -s amplicon_ref.fa ./outdir/amp_paired 10"<<endl<<endl;

        exit(0);
    }

    if(!fixed_seed){
	    rand_seed=(unsigned int) time(NULL); 
    }
    srand(rand_seed);

    bool is_pairend_read=false;
    char* seq_file= argv[i]; //"w:/vmLinux/vsProjects/art/art/release/test.fa";
    string out_file_prefix=argv[i+1];//"w:/vmLinux/vsProjects/art/art/release/test454";
    double x_fold  = atof(argv[i+2]);
    string num="";
    int fsize_mean=0;
    double fsize_std=0;

    if (k==6){
        num="1";
        is_pairend_read=true;
//        art::ini_read_pair_rand(abs(atoi(argv[i+3])),fabs(atof(argv[i+4])));

	fsize_mean=abs(atoi(argv[i+3]));
       	fsize_std=fabs(atof(argv[i+4]));
        art::ini_read_pair_rand(fsize_mean,fsize_std,rand_seed);

        if(art::gaussain_mean <= MINIMUM_FLAGMENT_SIZE){
            cerr<<"Error: the mean fragment length is shorter than the minimum fragment size defined ("<< MINIMUM_FLAGMENT_SIZE <<")"<<endl;
            exit(1);
        }
    }

    if(amplicon_pair){
        is_pairend_read=true;
//        art::ini_read_pair_rand(500,20,rand_seed);
    }

    string alnfasta=out_file_prefix+num+".aln";
    string fqfile=out_file_prefix+num+".fq";
    string statfile=out_file_prefix+".stat";

    string samfile=out_file_prefix+".sam";
    ofstream SAMFILE;
    if(sam_out) {
	    SAMFILE.open(samfile.c_str(),ios::binary); 
	    if(!SAMFILE.is_open()) { cerr<<"Can not open output file: "<<samfile<<endl; exit(0); }
    }

    ofstream FQFILE(fqfile.c_str(),ios::binary);
    if(!FQFILE.is_open()) { cout<<"can not open output file: "<<fqfile<<endl; exit(0); }

    ofstream ALNFILE;
    if(aln_out) {
	    ALNFILE.open(alnfasta.c_str(),ios::binary);
	    if(!ALNFILE.is_open()) { cout<<"can not open output file: "<<alnfasta<<endl; exit(0); }
	    ALNFILE<<"##ART_454"<<endl;
    }

    ofstream STATFILE(statfile.c_str(),ios::binary);
    if(!STATFILE.is_open()) { cout<<"can not open output file: "<<statfile<<endl; exit(0); }

    read_profile::set_err_prob();
    read_profile qdist; 
    if(profile_name.empty()){
	    qdist.default_profile(titanium); 
    }
    else{
	    qdist.user_profile(profile_name+"/qual_1st_profile", profile_name+"/qual_mc_profile", profile_name+"/indel_error_profile", profile_name+"/length_dist");
    }


    //caluate CPUT time
    clock_t start, end;
    double cpu_time_used;
    start = clock();

    samHeader sH;
    sH.getRefseqID(seq_file);
    sH.ID="02";
    sH.PN="ART_454"; 
    for(int j=0;j < argc; ++j) { sH.CL.append(argv[j]);sH.CL.append(" "); }
    if(aln_out){
	    sH.printAlnHeader(ALNFILE);
    }

    if(sam_out){
	    sH.printHeader(SAMFILE);
    }

    samRead sR;
    string srID;

    int read_len  = 0;
    vector<short> qual;
    readSeqFile seq_reader(seq_file);
    string id;
    art a_art; 
    seqRead a_read;
    string aln_read,aln_ref;
    ostringstream osID;
    int num_seq=0;
    string read_id;

    unsigned long total_read_count=0;

    string alnfasta2=out_file_prefix+"2.aln";
    string fqfile2=out_file_prefix+"2.fq";

    //paired-end reads simulaiton
    if(is_pairend_read){
	samRead sR2; 
	sR.rNext="=";
	sR2.rNext="=";


        ofstream FQFILE2(fqfile2.c_str(),ios::binary);
        if(!FQFILE2.is_open()) { cout<<"can not open output file: "<<fqfile2<<endl; exit(0); }

	ofstream ALNFILE2;
	if(aln_out) {
	       	ALNFILE2.open(alnfasta2.c_str(),ios::binary);
	       	if(!ALNFILE2.is_open()) { cout<<"can not open output file: "<<alnfasta2<<endl; exit(0); }
	       	ALNFILE2<<"##ART_454"<<endl;
	       	sH.printAlnHeader(ALNFILE2);
       	}

        seqRead a_read_2;
        vector<short> qual_2;
        string read_id_2;
        string aln_read_2,aln_ref_2;
       
        while(seq_reader.next_seq(id,a_art.ref_seq)){
	    std::replace( a_art.ref_seq.begin(), a_art.ref_seq.end(), 'U', 'T'); //replace U with T
            num_seq++;
            long t_num_base=(long) a_art.ref_seq.size();
            vector<int> depth_coverage(t_num_base,0);
            vector<int> start_coverage(t_num_base,0); //#read starting on the same position 

            long base_count=t_num_base; 
            int the_fold=x_fold;
            a_art.init_fast();
            while(the_fold>0){
                a_read.clear();

	    	int cc_try=0; //debug info
	       	if(amplicon){
		       	read_len=qdist.get_ran_read_len();
		       	a_art.ini_set(read_len); 
			if(!a_art.amp_pair_read(a_read, a_read_2)){
				cerr<<"Warning: skip the amplicon "<<id<<" as it is too short"<<endl; 
				break;
			}
	       	}
	       	else{
		       	do {
			       	read_len=qdist.get_ran_read_len();
			       	a_art.ini_set(read_len);
			       	//debug info
				if(cc_try>=100){
				       	cerr<<" was failed"<<endl;
			       	}
			       	cc_try++;
			       	if(cc_try>=100){ cerr<<"try "<<cc_try;}

			}while(!a_art.next_pair_read(a_read, a_read_2));
	       	}

		if(a_read.seq_ref.find('n',0)!=string::npos || a_read_2.seq_ref.find('n',0)!=string::npos) { 
			base_count-=read_len;
		       	continue;
	       	}

		//debug info
		if(cc_try>=100){ cerr<<" was successed"<<endl; }

                qual.clear();
                qual_2.clear();
                vector<string> aln, aln2;

		short start_cyc=0; //reset start_cyc to zero 
		bool ok=true;
                if(a_read.is_plus_strand){ 
                  qdist.get_read_qual_fast(a_art.homo_plus, a_read.bpos, a_read.seq_ref, a_read.seq_read, aln, qual, start_cyc, num_flow_cycles);
                  ok=qdist.get_read_qual_fast(a_art.homo_minus, a_read_2.bpos, a_read_2.seq_ref, a_read_2.seq_read, aln2, qual_2, start_cyc, num_flow_cycles);
                }
                else{
                  qdist.get_read_qual_fast(a_art.homo_minus, a_read.bpos, a_read.seq_ref, a_read.seq_read, aln, qual, start_cyc, num_flow_cycles);
                  ok=qdist.get_read_qual_fast(a_art.homo_plus, a_read_2.bpos, a_read_2.seq_ref, a_read_2.seq_read, aln2, qual_2, start_cyc, num_flow_cycles);
                }
		if(!ok) {

			if(debug){
			       	cerr <<"Info: "<<id<<", 1st read used all flowcycles, retry ..."<<endl;
			}
			continue;
		}

                //statistics
                start_coverage[a_read.bpos]+=1;
                start_coverage[a_read_2.bpos]+=1;
                for(int k=0; k<a_read.seq_ref.length(); k++){
                    depth_coverage[a_read.bpos+k]+=1;
                }
                for(int k=0; k<a_read_2.seq_ref.length(); k++){
                    depth_coverage[a_read_2.bpos+k]+=1;
                }

                //output
		total_read_count++;
                osID.str("");
                osID<<id<<'_'<<total_read_count; 
                read_id = osID.str();

		srID=read_id;
                read_id_2=read_id+"-2";
                read_id+="-1";
                //print first read

                FQFILE<<"@"<<read_id<<endl<<a_read.seq_read<<endl<<"+"<<endl;
                for(size_t k=0; k<qual.size(); k++){
                    FQFILE<<(char)(qual[k]+33);
                }
                FQFILE<<endl;

		if(aln_out){
		       	ALNFILE<<">"<<id<<"\t"<<read_id<<"\t"<<a_read.bpos;
		       	if(a_read.is_plus_strand) ALNFILE<<"\t+\n";
		       	else ALNFILE<<"\t-\n";
		       	ALNFILE<<aln[0]<<endl<<aln[1]<<endl;
		}	
                //print second read
                FQFILE2<<"@"<<read_id_2<<endl<<a_read_2.seq_read<<endl<<"+"<<endl;
                for(size_t k=0; k<qual_2.size(); k++){
                    FQFILE2<<(char)(qual_2[k]+33);
                }
                FQFILE2<<endl;
    
		if(aln_out) {
		       	ALNFILE2<<">"<<id<<"\t"<<read_id_2<<"\t"<<a_read_2.bpos;
		       	if(a_read_2.is_plus_strand) ALNFILE2<<"\t+\n";
		       	else ALNFILE2<<"\t-\n";
		       	ALNFILE2<<aln2[0]<<endl<<aln2[1]<<endl;
		}

		if(sam_out){
		       	sR.qname=srID;
		       	sR.rname=id;

		       	sR2.qname=srID;
		       	sR2.rname=id;

			sR.seq=a_read.seq_read;
		       	sR.qual.resize(qual.size());
		       	for(size_t k=0; k<qual.size(); k++){
			       	sR.qual[k]=(char)(qual[k]+33);
		       	}
		
			sR2.seq=a_read_2.seq_read;
		       	sR2.qual.resize(qual_2.size());
		       	for(size_t k=0; k<qual_2.size(); k++){
			       	sR2.qual[k]=(char)(qual_2[k]+33);
		       	}


			sR.flag=0x01 | 0x02 | 0x40;
		       	sR2.flag=0x01 | 0x02 | 0x80;

		       	if(a_read.is_plus_strand){
				       	sR.pos=a_read.bpos+1;
				       	sR.flag =sR.flag | 0x20;
				       	sR2.flag =sR2.flag |0x10;
				       	sR2.pos=a_art.ref_seq.size()-(a_read_2.bpos+a_read_2.seq_read.size()-1);
					sR2.reverse_comp();
		       	}
		       	else{
				       	sR.pos=a_art.ref_seq.size()-(a_read.bpos+a_read.seq_read.size()-1);
				       	sR.reverse_comp();
				       	sR2.pos=a_read_2.bpos+1;
				       	sR.flag = sR.flag | 0x10;
				       	sR2.flag = sR2.flag | 0x20;
		       	}

/*
		       	if(a_read.is_plus_strand){
			       	sR.pos=a_read.bpos+1;
			}
			else{
			       	sR.pos=a_art.ref_seq.size()-(a_read.bpos+a_read.seq_read.size()-1);
			       	sR.flag =sR.flag | 0x10;
			       	sR.reverse_comp();
			       	sR2.flag =sR2.flag |0x20;
			}

		       	if(a_read_2.is_plus_strand){
			       	sR2.pos=a_read_2.bpos+1;
			}
			else{
				sR2.pos=a_art.ref_seq.size()-(a_read_2.bpos+a_read_2.seq_read.size()-1);
			       	sR2.flag =sR2.flag |0x10;
			       	sR2.reverse_comp();
			       	sR.flag =sR.flag | 0x20;
			}
*/
		       	sR.getCigar(aln[0],aln[1],use_cigarM);
			sR2.getCigar(aln2[0],aln2[1],use_cigarM);

		       	sR.pNext=sR2.pos;
		       	sR2.pNext=sR.pos;
		       	if(sR2.pos>sR.pos){
			       	sR.tLen=sR2.pos+a_read_2.seq_read.size()-sR.pos;
			       	sR2.tLen=-sR.tLen;
		       	}
		       	else{
			       	sR2.tLen=sR.pos+a_read.seq_read.size()-sR2.pos;
			       	sR.tLen=-sR2.tLen;
		       	} 
			sR.printRead(SAMFILE);
			sR2.printRead(SAMFILE);
		}

		if(amplicon){
                    the_fold-=1;
		}
		else{
		       	base_count-=a_art.read_len;
		       	if(base_count<=0){
			       	the_fold-=1;
			       	base_count=t_num_base;
		       	}
		}
            }
            //output stat
            for (long k=0; k<t_num_base; k++){
                STATFILE<<k<<"\t"<<start_coverage[k]<<"\t"<<depth_coverage[k]<<endl;
            }
        }
        FQFILE2.close();
       	if(aln_out) { ALNFILE2.close();}
    }
    else{
	    //single-end reads simulaiton
        while(seq_reader.next_seq(id,a_art.ref_seq)){
	    std::replace( a_art.ref_seq.begin(), a_art.ref_seq.end(), 'U', 'T'); //replace U with T
            num_seq++;
            long t_num_base=(long) a_art.ref_seq.size();
            vector<int> depth_coverage(t_num_base,0);
            vector<int> start_coverage(t_num_base,0); //#read starting on the same position 
            long base_count=t_num_base; 
            int the_fold=x_fold;
            a_art.init_fast();

            while(the_fold>0){
                a_read.clear();
                read_len=qdist.get_ran_read_len();
                a_art.ini_set(read_len);
		if(amplicon) a_art.amp_read(a_read);
		else a_art.next_read(a_read);
		read_len=a_read.seq_ref.size();
                qual.clear();

                vector<string> aln;
		short start_cyc=0; //alway start with flow cycle zero for single-end reads

		if(a_read.seq_ref.find('n',0)!=string::npos) { 
			base_count-=read_len;
		       	continue;
	       	}
                
                if(a_read.is_plus_strand){ 
                  qdist.get_read_qual_fast(a_art.homo_plus, a_read.bpos, a_read.seq_ref, a_read.seq_read, aln, qual,start_cyc, num_flow_cycles);
                }
                else{
                  qdist.get_read_qual_fast(a_art.homo_minus, a_read.bpos, a_read.seq_ref, a_read.seq_read, aln, qual,start_cyc, num_flow_cycles);
                }

                //statistics
                start_coverage[a_read.bpos]+=1;
                for(int k=0; k<read_len; k++){
                    depth_coverage[a_read.bpos+k]+=1;
                }

                //output
                osID.str("");
                osID<<id<<'-'<<a_read.bpos<<'-'<<start_coverage[a_read.bpos]; 
                read_id = osID.str();

                FQFILE<<"@"<<read_id<<endl<<a_read.seq_read<<endl<<"+"<<endl;
                for(size_t k=0; k<qual.size(); k++){
                    FQFILE<<(char)(qual[k]+33);
                }
                FQFILE<<endl;
	
		if(aln_out){
		       	ALNFILE<<">"<<id<<"\t"<<read_id<<"\t"<<a_read.bpos;
		       	if(a_read.is_plus_strand) ALNFILE<<"\t+\n";
		       	else ALNFILE<<"\t-\n";
		       	ALNFILE<<aln[0]<<endl<<aln[1]<<endl;
		}

		if(sam_out){
		       	sR.qname=read_id;
		       	sR.rname=id;
		       	sR.flag =0;

			sR.seq=a_read.seq_read;
		       	sR.qual.resize(qual.size());
		       	for(size_t k=0; k<qual.size(); k++){
			       	sR.qual[k]=(char)(qual[k]+33);
		       	}

		       	if(a_read.is_plus_strand){
			       	sR.pos=a_read.bpos+1;
		       	}
		       	else{
			       	sR.flag = 0x10;
			       	sR.pos=a_art.ref_seq.size()-(a_read.bpos+a_read.seq_read.size()-1);
			       	sR.reverse_comp();
		       	}

		       	sR.getCigar(aln[0],aln[1],use_cigarM);
			sR.printRead(SAMFILE);
		}

		if(amplicon){
                    the_fold-=1;
		}
		else{
		       	base_count-=read_len;
		       	if(base_count<=0){
			       	the_fold-=1;
			       	base_count=t_num_base;
		       	}
	       	}
            }
            //output stat
            for (long k=0; k<t_num_base; k++){
                STATFILE<<k<<"\t"<<start_coverage[k]<<"\t"<<depth_coverage[k]<<endl;
            }
        }
    }

    FQFILE.close();
    if(aln_out){ ALNFILE.close(); }
    STATFILE.close();

    if(sam_out){ SAMFILE.close(); }

    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    
    if(!is_pairend_read){
	    if(amplicon){
		    cout << "            Amplicon 5'-end sequencing with single-end reads" << endl << endl;
	    }
	    else{
		    cout << "                          Single-end simulation" << endl << endl;
	    }
    } else {
	    if(amplicon){
		    cout << "           Amplicon two-end sequencing with paired-end reads" << endl << endl;
	    }
	    else{
		    cout << "                          Paired-end simulation" << endl << endl;
	    }
    }
    cout<<"Total CPU time used: "<<cpu_time_used<<endl<<endl;
    cout<<"The random seed for the run:   "<<rand_seed<<endl<<endl;

    cout << "Parameters Settings" << endl; 
    cout << "\tnumber of flow cycles:            " <<num_flow_cycles << endl;
    if(amplicon_pair){
	    cout << "\t# read pairs per amplion:         " << x_fold << endl;
    }
    else if (amplicon){
	    cout << "\t# reads per amplion:              " << x_fold << endl;
    }
    else{
	    cout << "\tfold of read coverage:            " << x_fold << "X" << endl;
    }
    if(is_pairend_read && !amplicon){
	    cout << "\tfragment length"<<endl;
	    cout << "\t\tmean:     " << fsize_mean << endl;
	    cout << "\t\tstd:      " << fsize_std << endl;
    }
    cout <<endl<<"454 Profile for Simulation" << endl; 
    if(profile_name.empty()){
	    if(titanium) cout << "\tthe built-in GS-FLX Titanium profile"<<endl; 
	    else cout << "\tthe built-in GS-FLX profile"<<endl; 
    }
    else{
	    cout << "\tthe profile provided: "<<profile_name<<endl; 
	    cout<<"\t\t"<<profile_name+"/qual_1st_profile"<<endl;
	    cout<<"\t\t"<<profile_name+"/qual_mc_profile"<<endl;
	    cout<<"\t\t"<<profile_name+"/indel_error_profile"<<endl;
	    cout<<"\t\t"<<profile_name+"/length_dist"<<endl;
    }
    cout <<endl<<"Output Files" << endl << endl;

    if(is_pairend_read){
	    cout << "  FASTQ Sequence Files:" << endl; 
	    cout << "\t the 1st reads: " << fqfile << endl;
	    cout << "\t the 2nd reads: " << fqfile2 << endl << endl;
	    if(aln_out){
		    cout << "  ALN Alignment Files:" << endl; 
		    cout << "\t the 1st reads: " << alnfasta << endl;
		    cout << "\t the 2nd reads: " << alnfasta2 << endl << endl;;
	    }
    } else {
	    cout << "  FASTQ Sequence File:" << endl; 
	    cout << "\t" << fqfile << endl << endl;
	    if (aln_out){
		    cout << "  ALN Alignment File:" << endl; 
		    cout << "\t" << alnfasta << endl << endl;
	    }
    }
    if(sam_out){
	    cout << "  SAM Alignment File:" << endl; 
	    cout << "\t" << samfile << endl << endl;
    }

	    cout << "  Read Coverage File:" << endl; 
	    cout << "\t" << statfile<< endl << endl;
}


