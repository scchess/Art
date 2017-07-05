/*
 * >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 * ART -- Artificial Read Transcription, ART_SOLiD
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
#include "readSeqFile.h"
#include "samRead.h"
#include "art.h"

using namespace std;
#define VERSION "1.3.2"
#define PRGNAME "art_SOLiD"
double SOLiDread::prob_err[max_qual_value];
//bool parse_arg(int num, char* arg);
//vector<double> SOLiDread::ins_rate;
//vector<double> SOLiDread::del_rate;
//vector<double> SOLiDread::sub_rate;
gsl_rng* art::gsl_R;
int art::gaussain_mean;
double art::gaussain_sigma;

int main(int argc, char* argv[]){
    cout << "================================================================="<<endl;
    cout << "                    ART_SOLiD (Version 1.3.3)                    "<<endl; 
    cout <<"         Simulation of Applied Biosystems' SOLiD Sequencing       "<<endl; 
    cout << "   Copyright (c) 2008-2015, Weichun Huang. All Rights Reserved.  "<<endl; 
    cout << "================================================================="<<endl<<endl;

    bool amplicon=false;
    char amp_read_type='s'; //type s:single-end, m:matepair, p: paired-end 
    bool is_pairend_read=false;
    bool is_mate_pair=false;
    bool use_cigarM = false;
    bool sam_out=false;

    long masked_reads_count=0;
    bool fixed_seed = false;
    unsigned int rand_seed = 0;

    string profile_name = "";
    double error_scale_factor=1;
    int i=1;
    for(;i<argc;++i){
	    char* pch = argv[i];
	    if( *pch != '-' || *(pch+1) == '\0') break;
	    while(*++pch){
		    switch(*pch){
/*
			    case 'v':
			    case 'V':
				    cout << "================================================================="<<endl;
				    cout << "                    ART_SOLiD (Version 1.3.2)                    "<<endl; 
				    cout <<"         Simulation of Applied Biosystems' SOLiD Sequencing       "<<endl; 
				    cout << "   Copyright (c) 2008-2014, Weichun Huang. All Rights Reserved.  "<<endl; 
				    cout << "================================================================="<<endl<<endl;
				    exit(1);
				    break;
*/
			    case 'p':
			    case 'P':
				    if(i<argc) profile_name=argv[++i]; 
				    else {
					    cerr<<"Error: no profile provided"<<endl;
					    exit(1);
				    }
				    break;
			    case 'A':
				    amplicon=true;
				    if(i<argc) amp_read_type=argv[++i][0];
				    else {
					    cerr<<"Error: no amplicon sequencing read type provided "<<endl;
					    exit(1);
				    } 
				    switch(amp_read_type){
					    case 's':
					    case 'm':
					    case 'p':
						    break;
					    default:
						    cerr<<"Error: undefined amplicon read type: "<<amp_read_type<<endl;
						    cerr<<"       Read type must be: s (single-end), m (matepair), or p (paired-end)"<<endl;
						    exit(1);

				    }
				    break;
			    case 'M':
				    use_cigarM=true;
				    break;
			    case 's':
			    case 'S':
				    sam_out=true;
				    break;
			    case 'f':
			    case 'F': 
				    if(i<argc) error_scale_factor=atof(argv[++i]);
				    else {
					    cerr<<"Error: no scaling factor provided"<<endl;
					    exit(1);
				    }
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

    int k=argc-i+1;
    if(k !=5 && k !=6 && k !=7 &&  k !=8){
        cout<<"===== USAGES ====="<<endl<<endl;
        cout<<"SINGLE-END (F3 READ) SIMULATION"<<endl;
        cout<<"\t"<< PRGNAME <<" [ options ] <INPUT_SEQ_FILE> <OUTPUT_FILE_PREFIX> <LEN_READ> <FOLD_COVERAGE>"<<endl<<endl;
//       	cout<<"     Example:"<<endl;
//       	cout<<"             "<<PRGNAME <<" -s seq_reference.fa ./outdir/single_dat 25 10"<<endl<<endl;
        cout<<"MATE-PAIR READS (F3-R3 PAIR) SIMULATION"<<endl;
        cout<<"\t"<< PRGNAME <<" [ options ] <INPUT_SEQ_FILE> <OUTPUT_FILE_PREFIX> <LEN_READ> <FOLD_COVERAGE> <MEAN_FRAG_LEN> <STD_DEV>"<<endl<<endl;
//        cout<<"     Example:"<<endl;
//       	cout<<"             1) simulation of matepair reads"<<endl;
//       	cout<<"             "<<PRGNAME <<" -s seq_reference.fa ./outdir/matepair_dat 32 10 1000 50"<<endl<<endl;
//       	cout<<"             2) simulation of matepair reads with a fixed random seed"<<endl;
//       	cout<<"             "<<PRGNAME <<" -r 777 -s seq_reference.fa ./outdir/matepair_fs 50 10 1000 50"<<endl<<endl;
        cout<<"PAIRED-END READS (F3-F5 PAIR) SIMULATION"<<endl;
        cout<<"\t"<< PRGNAME <<" [ options ] <INPUT_SEQ_FILE> <OUTPUT_FILE_PREFIX> <LEN_READ_F3> <LEN_READ_F5> <FOLD_COVERAGE> <MEAN_FRAG_LEN> <STD_DEV>"<<endl<<endl;
//       cout<<"     Example:"<<endl;
//       	cout<<"             "<<PRGNAME <<" -s seq_reference.fa ./outdir/paired_dat 33 25 10 250 10"<<endl<<endl;
        cout<<"AMPLICON SEQUENCING SIMULATION"<<endl;
        cout<<"\t"<< PRGNAME <<" [ options ] -A s <INPUT_SEQ_FILE> <OUTPUT_FILE_PREFIX> <LEN_READ> <READS_PER_AMPLICON>"<<endl;
        cout<<"\t"<< PRGNAME <<" [ options ] -A m <INPUT_SEQ_FILE> <OUTPUT_FILE_PREFIX> <LEN_READ> <READ_PAIRS_PER_AMPLICON>"<<endl;
        cout<<"\t"<< PRGNAME <<" [ options ] -A p <INPUT_SEQ_FILE> <OUTPUT_FILE_PREFIX> <LEN_READ_F3> <LEN_READ_F5> <READ_PAIRS_PER_AMPLICON>"<<endl<<endl;
//       	cout<<"     Examples:"<<endl;
//       	cout<<"             1) amplicon sequencing with single-end reads"<<endl;
//       	cout<<"             "<<PRGNAME <<" -A s -s amp_reference.fa ./outdir/amp_single_end_dat 25 100"<<endl<<endl;
//       	cout<<"             2) amplicon sequencing with matepair reads"<<endl;
//       	cout<<"             "<<PRGNAME <<" -A m -s amp_reference.fa ./outdir/amp_single_end_dat 50 50"<<endl<<endl;
//       	cout<<"             3) amplicon sequencing with paired-end reads"<<endl;
//       	cout<<"             "<<PRGNAME <<" -A p -s amp_reference.fa ./outdir/amp_single_end_dat 35 25 50"<<endl<<endl;
        cout<<"===== PARAMETERS ====="<<endl<<endl;
       	cout<<"INPUT_SEQ_FILE            -  filename of DNA/RNA reference sequences in FASTA format"<<endl;
       	cout<<"OUTPUT_FILE_PREFIX        -  prefix or directory for all output read data files"<<endl;
       	cout<<"FOLD_COVERAGE             -  fold of read coverage over the reference sequences "<<endl;
       	cout<<"LEN_READ                  -  length of F3/R3 reads"<<endl;
       	cout<<"LEN_READ_F3               -  length of F3 reads for paired-end read simulation"<<endl;
       	cout<<"LEN_READ_F5               -  length of F5 reads for paired-end read simulation"<<endl;
       	cout<<"READS_PER_AMPLICON        -  number of reads per amplicon"<<endl;
       	cout<<"READ_PAIRS_PER_AMPLICON   -  number of read pairs per amplicon"<<endl;
       	cout<<"MEAN_FRAG_LEN             -  mean DNA/RNA fragment size for matepair/paired-end read simulation"<<endl; 
	cout<<"STD_DEV                   -  standard deviation of the DNA/RNA fragment sizes for matepair/paired-end read simulation"<<endl<<endl; 

        cout<<"===== OPTIONAL PARAMETERS ====="<< endl<<endl;
        cout<<"-A specify the read type for amplicon sequencing simulation (s:single-end, m: matepair, p: paired-end)"<<endl;
        cout<<"-M indicate to use CIGAR 'M' instead of '=/X' for alignment match/mismatch"<<endl;
        cout<<"-s indicate to generate a SAM alignment file"<<endl;
        cout<<"-r specify the random seed for the simulation"<<endl;
        cout<<"-f specify the scale factor adjusting error rate (e.g., -f 0 for zero-error rate simulation)"<<endl;
        cout<<"-p specify user's own read profile for simulation"<<endl<<endl;
        cout<<"===== EXAMPLES ====="<< endl<<endl;
        cout<<" 1) singl-end 25bp reads simulation at 10X coverage"<<endl;
       	cout<<"\t"<<PRGNAME <<" -s seq_reference.fa ./outdir/single_dat 25 10"<<endl;
        cout<<" 2) singl-end 75bp reads simulation at 20X coverage with user's error profile"<<endl;
       	cout<<"\t"<<PRGNAME <<" -s -p ../SOLiD_profiles/profile_pseudo ./seq_reference.fa ./dat_userProfile 75 20"<<endl;
       	cout<<" 3) matepair 35bp (F3-R3) reads simulation at 20X coverage with DNA/RNA MEAN fragment size 2000bp and STD 50"<<endl;
       	cout<<"\t"<<PRGNAME <<" -s seq_reference.fa ./outdir/matepair_dat 35 20 2000 50"<<endl;
       	cout<<" 4) matepair reads simulation with a fixed random seed"<<endl;
       	cout<<"\t"<<PRGNAME <<" -r 777 -s seq_reference.fa ./outdir/matepair_fs 50 10 1500 50"<<endl;
        cout<<" 5) paired-end reads (75bp F3, 35bp F5) simulation with the MEAN fragment size 250 and STD 10 at 20X coverage"<<endl;
       	cout<<"\t"<<PRGNAME <<" -s seq_reference.fa ./outdir/paired_dat 75 35 50 250 10"<<endl;
       	cout<<" 6) amplicon sequencing with 25bp single-end reads at 100 reads per amplicon"<<endl;
       	cout<<"\t"<<PRGNAME <<" -A s -s amp_reference.fa ./outdir/amp_single 25 100"<<endl;
       	cout<<" 7) amplicon sequencing with 50bp matepair reads at 80 read pairs per amplicon"<<endl;
       	cout<<"\t"<<PRGNAME <<" -A m -s amp_reference.fa ./outdir/amp_matepair 50 80"<<endl;
       	cout<<" 8) amplicon sequencing with paired-end reads (35bp F3, 25bp F5 reads) at 50 pairs per amplicon"<<endl;
       	cout<<"\t"<<PRGNAME <<" -A p -s amp_reference.fa ./outdir/amp_pair 35 25 50"<<endl<<endl;
//        cout<<"-v print out version information"<<endl;

        exit(1);
    }

    if(!fixed_seed){
	    rand_seed=(unsigned int) time(NULL); 
    }
    srand(rand_seed);

    bool mask_n=true; 
    short max_num_n=1; 
    int len_ref_id=250;

    //caluate CPUT time
     clock_t start, end;
     double cpu_time_used;
     start = clock();

    char* seq_file= argv[i];
    string out_file_prefix=argv[i+1];
    int read_len  = atoi(argv[i+2]);
    if(read_len<=0){ cerr<<"Error: read length must be > 0"<<endl; exit(1); }
    int read_len_R3F5  =0;
    double x_fold  =0;
    int fsize_mean=0;
    double fsize_std=0;
    string num="";

    if (k==7 || k==8){
        is_pairend_read=true;
	if(k==7){
	       	num="_R3";
	       	is_mate_pair=true;
	       	read_len_R3F5 =read_len;
	       	x_fold  = atof(argv[i+3]);
	       	fsize_mean=abs(atoi(argv[i+4]));
	       	fsize_std=fabs(atof(argv[i+5]));
       	}
	else{
	       	num="_F5";
	       	is_mate_pair=false;
	       	read_len_R3F5 =atoi(argv[i+3]);
		if(read_len_R3F5<=0){
			cerr<<"Error: read lenght must be > 0"<<endl; exit(1);
		}
	       	x_fold  = atof(argv[i+4]);
	       	fsize_mean=abs(atoi(argv[i+5]));
	       	fsize_std=fabs(atof(argv[i+6]));
       	}
//        if(fixed_seed){
	       	art::ini_read_pair_rand(fsize_mean,fsize_std, rand_seed);
//	}
//	else{
//		art::ini_read_pair_rand(abs(atoi(argv[i+4])),fabs(atof(argv[i+5])));
//	}
        if(art::gaussain_mean<=read_len){
            cerr<<"Error: the read length must be shorter than the mean flagment length specified"<<endl;
            exit(1);
        }
    }
    else{
        is_pairend_read=false;
        is_mate_pair=false;
	if(k==6){
	       	read_len_R3F5 =atoi(argv[i+3]);
	       	x_fold  = atof(argv[i+4]);
	       	if(!amplicon){
		       	cerr<<"Error: wrong usage"<<endl;
			exit(1);
		}
	}
	else{
	       	x_fold  = atof(argv[i+3]);
	}
       	if(amplicon){
	       	if(amp_read_type=='m'){
		       	is_pairend_read=true;
		       	is_mate_pair=true;
		       	num="_R3";
		       	read_len_R3F5 =read_len;
	       	}
	       	else if(amp_read_type=='p'){
		       	num="_F5";
		       	is_pairend_read=true;
	       	}
       	}
    }

    if(read_len >35 || read_len_R3F5 >35){
	    if(profile_name.empty()){
		    if(read_len<=75 && read_len_R3F5 <=75){
			    profile_name="pseudo";
			    cerr<<"Warning: you are using the 75bp error profile for testing only."<<endl<<endl;  
		    }
		    else{
			    cerr<<"Error: the read length exceeds the max length (75bp) with the built-in error profile"<<endl;
			    exit(1);
		    }
	    }
    }

    string seqfasta=out_file_prefix+num+".fa";
    string qualfasta=out_file_prefix+num+".qual";
    string alnfasta=out_file_prefix+num+".map";
    string fqfile=out_file_prefix+num+".fq";

    string samfile=out_file_prefix+".sam";
    ofstream SAMFILE;
    if(sam_out) {
	    SAMFILE.open(samfile.c_str(),ios::binary); 
	    if(!SAMFILE.is_open()) { cerr<<"Can not open output file: "<<samfile<<endl; exit(0); }
    }

    ofstream FQFILE(fqfile.c_str(),ios::binary);
    if(!FQFILE.is_open()) { cout<<"can not open output file: "<<fqfile<<endl; exit(0); }

    ofstream ALNFILE(alnfasta.c_str(),ios::binary);
    if(!ALNFILE.is_open()) { cout<<"can not open output file: "<<alnfasta<<endl; exit(0); }
    ALNFILE<<"##ART_SOLiD\tread_length\t"<<read_len<<endl;

    samHeader sH;
    sH.getRefseqID(seq_file);
    sH.ID="03";
    sH.PN="ART_SOLiD"; 
    for(int i=0;i < argc; ++i) { sH.CL.append(argv[i]);sH.CL.append(" "); }
    sH.printAlnHeader(ALNFILE);

    if(sam_out){
	    sH.printHeader(SAMFILE);
    }

    samRead sR;
    string srID;

    SOLiDread::set_err_prob();

    vector<short> qual;
    readSeqFile seq_reader(seq_file);
    string id;
    art a_art; 
    SOLiDread a_read(profile_name);
    a_read.ini_ran_qual(rand_seed);

    if(read_len>a_read.error_profile.size()){
        cerr<<"Error: the read length "<<read_len<<" exceeds the max length "<<a_read.error_profile.size()<<endl;
//        cerr<<"Error: the read length "<<read_len<<" exceeds the max length "<<qdist.qual_dist_first.size()<<endl<<endl;
        exit(1);
    }

    for(size_t i=0; i<a_read.cal_err_rate_1st.size(); i++){
      a_read.cal_err_rate_1st[i]*=error_scale_factor;
    }
    for(size_t i=0; i<a_read.cal_err_rate_2nd.size(); i++){
      a_read.cal_err_rate_2nd[i]*=error_scale_factor;
    }

/*
    a_read.set_rate(read_len,0.0001,2,a_read.ins_rate);
    a_read.set_rate(read_len,0.0001,2,a_read.del_rate);
    a_read.set_rate(read_len,0.028,2,a_read.sub_rate);
*/
    string aln_read,aln_ref;
    ostringstream osID;
    int num_seq=0;
    string read_id;
   
    string seqfasta2=out_file_prefix+"_F3.fa";
    string qualfasta2=out_file_prefix+"_F3.qual";
    string alnfasta2=out_file_prefix+"_F3.map";
    string fqfile2=out_file_prefix+"_F3.fq";
   
    unsigned long cc_num_read=1;
    if(is_pairend_read){
	samRead sR2; 
	sR.rNext="=";
	sR2.rNext="=";
//        ofstream SEQFILE2(seqfasta2.c_str(),ios::binary);
//        if(!SEQFILE2.is_open()) { cout<<"can not open output file: "<<seqfasta2<<endl; exit(0); }

//        ofstream QUALFILE2(qualfasta2.c_str(),ios::binary);
//        if(!QUALFILE2.is_open()) { cout<<"can not open output file: "<<qualfasta2<<endl; exit(0); }

        ofstream FQFILE2(fqfile2.c_str(),ios::binary);
        if(!FQFILE2.is_open()) { cout<<"can not open output file: "<<fqfile2<<endl; exit(0); }

        ofstream ALNFILE2(alnfasta2.c_str(),ios::binary);
        if(!ALNFILE2.is_open()) { cout<<"can not open output file: "<<alnfasta2<<endl; exit(0); }
       	ALNFILE2<<"##ART_SOLiD\tread_length\t"<<read_len<<endl;
       	sH.printAlnHeader(ALNFILE2);

        SOLiDread a_read_2(profile_name);
        a_read_2.ini_ran_qual(rand_seed);

        for(size_t i=0; i<a_read_2.cal_err_rate_1st.size(); i++){
          a_read_2.cal_err_rate_1st[i]*=error_scale_factor;
        }
        for(size_t i=0; i<a_read_2.cal_err_rate_2nd.size(); i++){
          a_read_2.cal_err_rate_2nd[i]*=error_scale_factor;
        }
/*
   	a_read_2.set_rate(read_len,0.0001,2,a_read.ins_rate);
        a_read_2.set_rate(read_len,0.0001,2,a_read.del_rate);
        a_read_2.set_rate(read_len,0.036,2,a_read.sub_rate);
*/
        vector<short> qual_2;
        string read_id_2;
        string aln_read_2,aln_ref_2;
        while(seq_reader.next_seq(id,a_art.ref_seq)){ 
//            size_t p1=id.find_first_of(' '); if(p1==string::npos) p1=10; size_t p2=id.find_first_of('\t'); if(p2==string::npos) p2=10;            p1=p1<p2?p1:p2; id=id.substr(0,p1); 
	    std::replace(a_art.ref_seq.begin(), a_art.ref_seq.end(), 'U', 'T'); //replace U with T

            istringstream isID; isID.str(id); isID>>id; id=id.substr(0,len_ref_id); 
            num_seq++;
            a_art.ini_set(read_len,read_len_R3F5);
            if(mask_n){ 
              a_art.mask_n_region(max_num_n);
            }
            long t_num_read=(long) a_art.ref_seq.size()/read_len*x_fold;
            while(t_num_read>0){
//generate SOLiD-like id
		int num_3rd=cc_num_read / 1000000 + 1;
		unsigned int num_3rd_rem=cc_num_read % 1000000;
		int num_1st=num_3rd_rem % 1000;
		int num_2nd=num_3rd_rem / 1000 + 1;
                osID<<num_3rd<<'_'<<num_2nd<<'_'<<num_1st;
                read_id = osID.str();
                osID.str("");

                a_read.clear();
                a_read_2.clear();
                //a_art.next_pair_read_indel(a_read, a_read_2); //need SOLiD profile with indel error rate 
		
		if(amplicon){
		       	if(is_mate_pair){
			       	a_art.amp_mate_read(a_read, a_read_2); 
			}
		       	else{
			       	a_art.amp_PE_read(a_read, a_read_2); 
			}
		}
		else{
		       	if(is_mate_pair){
			       	a_art.next_pair_read(a_read, a_read_2); 
			}
		       	else{
			       	a_art.next_PE_read(a_read, a_read_2); 
			}
		}
                if(mask_n){ 
                  if(a_read.is_plus_strand){
                    if(a_art.masked_pos.count(a_read.bpos)>0){
                      t_num_read-=2;
		      masked_reads_count+=1;
                      continue;
                    }
                  }
                  else{
                    size_t bpos1=a_art.ref_seq.size()-a_read.bpos-read_len_R3F5;
                    if(a_art.masked_pos.count(bpos1)>0){
		      masked_reads_count+=1;
                      t_num_read-=2;
                      continue;
                    }
                  }

                  if(a_read_2.is_plus_strand){
                    if(a_art.masked_pos.count(a_read_2.bpos)>0){
                      t_num_read-=2;
		      masked_reads_count+=1;
                      continue;
                    }
                  }
                  else{
                    size_t bpos2=a_art.ref_seq.size()-a_read_2.bpos-read_len;
                    if(a_art.masked_pos.count(bpos2)>0){
		      masked_reads_count+=1;
                      t_num_read-=2;
                      continue;
                    }
                  }
                }

		string cs_seq_1st, cs_seq_2nd;
		map<int,char> error_pos_1st, error_pos_2nd;
		map<int,char>::iterator it; 
		//convert base-space to color space, and incorporate sequencing errors
                //qual.clear();
                //qual_2.clear();
	    	a_read.convert_seq2cs(cs_seq_1st, qual, error_pos_1st, false) ; //the 1st is R3 
	    	a_read_2.convert_seq2cs(cs_seq_2nd, qual_2, error_pos_2nd) ;  //the 2nd is F3 
//                read_id_2=read_id+"-2";
//                read_id+="-1";
		srID=read_id;
                read_id_2=read_id+"_F3";
	       	read_id+=num;
//print first read
//                SEQFILE<<">"<<read_id<<endl<<a_read.seq_read<<endl; //<<a_read.seq_ref<<endl;
//                QUALFILE<<">"<<read_id<<endl;
//                copy(qual.begin(),qual.end(), ostream_iterator<short>(QUALFILE,"\t"));
//                QUALFILE<<endl;

                FQFILE<<"@"<<read_id<<endl<<'G'<<cs_seq_1st<<endl<<"+"<<endl;
                for(size_t k=0; k<qual.size(); k++){
                    FQFILE<<(char)(qual[k]+33);
                }
                FQFILE<<endl;

                ALNFILE<<id<<"\t"<<read_id<<"\t"<<a_read.bpos;
                if(a_read.is_plus_strand) ALNFILE<<"\t+";
                else ALNFILE<<"\t-";
		ALNFILE<<"\t"<<error_pos_1st.size();
		for(it=error_pos_1st.begin(); it!=error_pos_1st.end(); it++){
			ALNFILE<<"\t"<<it->first<<"\t"<<it->second<<cs_seq_1st[it->first];
		}
		ALNFILE<<endl;

		/*
                if(a_read.get_aln(aln_read,aln_ref)){
                    ALNFILE<<aln_ref<<endl<<aln_read<<endl;
                }
                else{
                    ALNFILE<<a_read.seq_ref<<endl<<a_read.seq_read<<endl;
                }
		*/
//print second read
//                SEQFILE2<<">"<<read_id_2<<endl<<a_read_2.seq_read<<endl; //<<a_read.seq_ref<<endl;
//                QUALFILE2<<">"<<read_id_2<<endl;
//                copy(qual_2.begin(),qual_2.end(), ostream_iterator<short>(QUALFILE2,"\t"));
//                QUALFILE2<<endl;

                FQFILE2<<"@"<<read_id_2<<endl<<'T'<<cs_seq_2nd<<endl<<"+"<<endl;
                for(size_t k=0; k<qual_2.size(); k++){
                    FQFILE2<<(char)(qual_2[k]+33);
                }
                FQFILE2<<endl;

                ALNFILE2<<id<<"\t"<<read_id_2<<"\t"<<a_read_2.bpos;
                if(a_read_2.is_plus_strand) ALNFILE2<<"\t+";
                else ALNFILE2<<"\t-";
		ALNFILE2<<"\t"<<error_pos_2nd.size();
		for(it=error_pos_2nd.begin(); it!=error_pos_2nd.end(); it++){
			ALNFILE2<<"\t"<<it->first<<"\t"<<it->second<<cs_seq_2nd[it->first];
		}
		ALNFILE2<<endl;
/*
                if(a_read_2.get_aln(aln_read_2,aln_ref_2)){
                    ALNFILE2<<aln_ref_2<<endl<<aln_read_2<<endl;
                }
                else{
                    ALNFILE2<<a_read_2.seq_ref<<endl<<a_read_2.seq_read<<endl;
                }
*/
                t_num_read-=2;
                cc_num_read+=1;

		if(sam_out){
		       	sR.qname=srID;
		       	sR.rname=id;

		       	sR2.qname=srID;
		       	sR2.rname=id;
//			sR.seq=a_read.seq_read;
		       	a_read.convert_cs2seq(cs_seq_1st, aln_read, 'g');
			sR.seq=aln_read;
		       	sR.qual.resize(qual.size());
		       	for(size_t k=0; k<qual.size(); k++){
			       	sR.qual[k]=(char)(qual[k]+33);
		       	}
		
//			sR2.seq=a_read_2.seq_read;
		       	a_read_2.convert_cs2seq(cs_seq_2nd, aln_read_2, 't');
			sR2.seq=aln_read_2;
		       	sR2.qual.resize(qual_2.size());
		       	for(size_t k=0; k<qual_2.size(); k++){
			       	sR2.qual[k]=(char)(qual_2[k]+33);
		       	}
			
			aln_ref=a_read.seq_ref;
			aln_ref_2=a_read_2.seq_ref;

			sR.flag=0x01 | 0x02 | 0x40;
		       	sR2.flag=0x01 | 0x02 | 0x80;

			if(is_mate_pair){
			       	if(a_read.is_plus_strand){
				       	sR.pos=a_read.bpos+1;
			       	}
			       	else{
				       	sR.pos=a_art.ref_seq.size()-(a_read.bpos+read_len-1);
				       	sR.flag =sR.flag | 0x10;
				       	sR.reverse_comp();
				       	sR2.flag =sR2.flag |0x20;
			       	} 
				if(a_read_2.is_plus_strand){
				       	sR2.pos=a_read_2.bpos+1;
			       	}
			       	else{
				       	sR2.pos=a_art.ref_seq.size()-(a_read_2.bpos+read_len-1);
				       	sR2.flag =sR2.flag |0x10;
				       	sR2.reverse_comp();
				       	sR.flag =sR.flag | 0x20;
			       	}
			}
			else{
			       	if(a_read.is_plus_strand){
				       	sR.pos=a_read.bpos+1;
				       	sR2.pos=a_art.ref_seq.size()-(a_read_2.bpos+read_len-1);
				       	sR2.flag =sR2.flag |0x10;
				       	sR2.reverse_comp();
				       	sR.flag =sR.flag | 0x20;
			       	}
			       	else{
				       	sR2.pos=a_read_2.bpos+1;
				       	sR.pos=a_art.ref_seq.size()-(a_read.bpos+read_len_R3F5-1);
				       	sR.flag =sR.flag | 0x10;
				       	sR.reverse_comp();
				       	sR2.flag =sR2.flag |0x20;
			       	} 
			}

		       	sR.getCigar(aln_ref,aln_read,use_cigarM);
			sR2.getCigar(aln_ref_2,aln_read_2,use_cigarM);

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

            }
        }
        FQFILE2.close();
        ALNFILE2.close();
    }
    else{
        while(seq_reader.next_seq(id,a_art.ref_seq)){
	    std::replace(a_art.ref_seq.begin(), a_art.ref_seq.end(), 'U', 'T'); //replace U with T
            istringstream isID; isID.str(id); isID>>id; id=id.substr(0,len_ref_id); 
            num_seq++;
            a_art.ini_set(read_len);
            if(mask_n){ 
              a_art.mask_n_region(max_num_n);
            }
            long t_num_read=(long) a_art.ref_seq.size()/read_len*x_fold;
            while(t_num_read>0){
//              osID<<num_seq<<fixed<<setfill('0')<< setw(10)<< t_num_read;
/*
                osID<<id<<'-'<<t_num_read;
                read_id = osID.str();
                osID.str("");
*/

//generate SOLiD-like id
		int num_3rd=cc_num_read / 1000000 + 1;
		unsigned int num_3rd_rem=cc_num_read % 1000000;
		int num_1st=num_3rd_rem % 1000;
		int num_2nd=num_3rd_rem / 1000 + 1;

                osID<<num_3rd<<'_'<<num_2nd<<'_'<<num_1st;
                read_id = osID.str();
                osID.str("");
		read_id +="_F3";

                a_read.clear();
                //a_art.next_read(a_read);
		if(amplicon) a_art.amp_read(a_read);
		else a_art.next_read(a_read);

                if(mask_n){ 
                  if(a_read.is_plus_strand){
                    if(a_art.masked_pos.count(a_read.bpos)>0){
                      t_num_read-=1;
		      masked_reads_count+=1;
                      continue;
                    }
                  }
                  else{
                    size_t bpos=a_art.ref_seq.size()-a_read.bpos-read_len;
                    if(a_art.masked_pos.count(bpos)>0){
                      t_num_read-=1;
		      masked_reads_count+=1;
                      continue;
                    }
                  }
                }

 		string cs_seq_1st;
		map<int,char> error_pos_1st;
		map<int,char>::iterator it; 
                //qua.clear();
	    	a_read.convert_seq2cs(cs_seq_1st, qual, error_pos_1st) ;
                FQFILE<<"@"<<read_id<<endl<<'T'<<cs_seq_1st<<endl<<"+"<<endl;
                for(size_t k=0; k<qual.size(); k++){
                    FQFILE<<(char)(qual[k]+33);
                }
                FQFILE<<endl;

                ALNFILE<<id<<"\t"<<read_id<<"\t"<<a_read.bpos;
                if(a_read.is_plus_strand) ALNFILE<<"\t+";
                else ALNFILE<<"\t-";
		ALNFILE<<"\t"<<error_pos_1st.size();
		for(it=error_pos_1st.begin(); it!=error_pos_1st.end(); it++){
			ALNFILE<<"\t"<<it->first<<"\t"<<it->second<<cs_seq_1st[it->first];
		}
		ALNFILE<<endl;

                t_num_read--;
                cc_num_read+=1;

		if(sam_out){
		       	sR.qname=read_id;
		       	sR.rname=id;

		       	a_read.convert_cs2seq(cs_seq_1st, aln_read, 't'); 
			aln_ref=a_read.seq_ref;
			sR.seq=aln_read;
		       	sR.qual.resize(qual.size());
		       	for(size_t k=0; k<qual.size(); k++){
			       	sR.qual[k]=(char)(qual[k]+33);
		       	}
		
			sR.flag=0;

		       	if(a_read.is_plus_strand){
			       	sR.pos=a_read.bpos+1;
			}
			else{
			       	sR.pos=a_art.ref_seq.size()-(a_read.bpos+read_len-1);
			       	sR.flag =0x10;
			       	sR.reverse_comp();
			}

		       	sR.getCigar(aln_ref,aln_read,use_cigarM);
			sR.printRead(SAMFILE);
		}
            }
        }
    }

    FQFILE.close();
    ALNFILE.close();
    if(sam_out){ SAMFILE.close(); }
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;

    if(amplicon){
	    if(!is_pairend_read){
		    cout << "          Amplicon single-end sequencing simulation" << endl << endl;
	    }
	    else if(is_mate_pair){
		    cout << "          Amplicon matepair sequencing simulation" << endl << endl;
	    }
	    else{
		    cout << "          Amplicon paired-end sequencing simulation" << endl << endl;
	    }
    }
    else{
	    if(!is_pairend_read){
		    cout << "                    Single-end simulation" << endl << endl;
	    }
	    else if(is_mate_pair){
		    cout << "                    Mate-Pair (F3-R3) simulation" << endl << endl;
	    }
	    else{
		    cout << "                    Paired-end (F3-F5) simulation" << endl << endl;
	    }
    }
    cout<<"Total CPU time used: "<<cpu_time_used<<endl<<endl;
    cout<<"The random seed for the run:   "<<rand_seed<<endl<<endl;

    cout << "Parameters Settings" << endl; 
    if (amplicon){
	    if(is_pairend_read)
		    cout << "\t# read pairs per amplion:\t" << x_fold << endl;
	    else
		    cout << "\t# reads per amplion:\t" << x_fold << endl;
    }
    else{
	    cout << "\tfold of read coverage:\t" << x_fold << "X" << endl;
    }
    if(is_pairend_read){
	    cout << "\tF3 read length:\t"<<read_len<<endl;
	    if(is_mate_pair) cout << "\tR3 read length:\t"<<read_len_R3F5<<endl;
	    else cout << "\tF5 read length:\t"<<read_len_R3F5<<endl;
	    if(!amplicon){
		    cout << "\tfragment length"<<endl;
		    cout << "\t\tmean:\t" << fsize_mean << endl;
		    cout << "\t\tstd:\t" << fsize_std << endl;
	    }
    }
    else{
	    cout << "\tread length:\t"<<read_len<<endl;
    }
    cout <<endl<<"SOLiD Error Profile for Simulation" << endl; 
    if(profile_name.empty()){
	    cout << "\tthe built-in 35bp error profile"<<endl; 
    }
    else if(profile_name=="pseudo"){
	    cout << "\tthe 75bp error profile for testing only"<<endl; 
    }
    else{
	    cout << "\tthe profile provided: "<<profile_name<<endl; 
    }
    cout <<endl<<"Output Files" << endl << endl;

    if(is_pairend_read){
	    cout << "  FASTQ Sequence Files:" << endl; 
	    cout << "\t the 1st reads: " << fqfile << endl;
	    cout << "\t the 2nd reads: " << fqfile2 << endl << endl;
//	    if(aln_out){
		    cout << "  MAP Alignment Files:" << endl; 
		    cout << "\t the 1st reads: " << alnfasta << endl;
		    cout << "\t the 2nd reads: " << alnfasta2 << endl << endl;;
//	    }
    } else {
	    cout << "  FASTQ Sequence File:" << endl; 
	    cout << "\t" << fqfile << endl << endl;
//	    if (aln_out){
		    cout << "  MAP Alignment File:" << endl; 
		    cout << "\t" << alnfasta << endl << endl;
//	    }
    }
    if(sam_out){
	    cout << "  SAM Alignment File:" << endl; 
	    cout << "\t" << samfile << endl << endl;
    }

    if(masked_reads_count){
	    cout << "NOTE: all genomic regions with 'N' were masked" << endl; 
	    if(is_pairend_read){
		    cout << "# discarded pairs of reads mapped to the masked regions:\t" <<masked_reads_count<<endl<<endl; 
	    }
	    else{
		    cout << "# discarded reads mapped to the masked regions:\t" <<masked_reads_count<<endl<<endl; 
	    }
    }

}

    //      ofstream outGlobal(outFile.c_str());
    //      if(!outGlobal.is_open()) { cout<<"can not open output file: "<<outFile<<endl; exit(0); }
    //      ostream_iterator <char, char, char_traits <char> > os(outGlobal,"");		

//
//bool parse_arg(int num, char* arg){
//    bool success=true;
//    int i=1;
//    for(;i<ARGC;++i){
//       	char* pch = ARGV[i];
//       	if( *pch != '-' || *(pch+1) == '\0') break;
//       	while(*++pch && success){
//            switch(*pch){
//                // Version Information
//            case 'v':
//            case 'V':
//                prtVersion = true;
//                break;
//                // Help Information
//            case 'h':
//            case 'H':
//                prtVersion = true;
//                prtUsage = true;
//                break;
//                // Unbuffered Output
//            case 't':
//            case 'T':
//                html_out=false;
//                break;
//            case 'b':
//            case 'B':
//                if(i<ARGC) browser=ARGV[++i]; 
//                else success=PrintErr("Invalid value for option: \"%c\"", *pch);
//                break;
//            case 'o':
//            case 'O':
//                if(i<ARGC) outFile=ARGV[++i]; 
//                else success=PrintErr("Invalid value for option: \"%c\"", *pch);
//                break;
//            case 'c':
//            case 'C':
//                if(i<ARGC) cfgFile=ARGV[++i]; 
//                else success=PrintErr("Invalid value for option: \"%c\"", *pch);
//                break;
//                // Error Reporting
//            default:
//                success=PrintErr("Unreconized switch, \"%c\"", *pch);
//                prtUsage = true;
//                prtVersion = true;
//                break;
//            }
//       	}
//    }
//
//    // Print Version and/or Usage Information and exit
//    if(prtVersion || prtUsage){
//       	if(prtUsage) printUsage();
//       	if(prtVersion) printVer();
//       	return false;
//    }
//
//    if(i>ARGC){ printUsage(); return false;}
//
//    aceFile=ARGV[i];
//    return success;;
//}

