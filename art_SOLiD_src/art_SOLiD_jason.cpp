/*
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
ART_SOLiD -- Assembly Read Transcriber
Copyright(c) 2008-2011 Weichun Huang, All Rights Reserved.
<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
*/
#include <iostream>
#include <sstream>
#include <string>
#include <time.h>
#include <algorithm>
#include <iomanip>
#include <ctime>
#include "readSeqFile.h"
#include "empdist.h"
#include "art.h"

using namespace std;
#define MAJOR_VERSION "0.5"
#define MINOR_VERSION "0"
//bool parse_arg(int num, char* arg);
double empdist::prob_err[HIGHEST_QUAL];
//vector<double> seqRead::ins_rate;
//vector<double> seqRead::del_rate;
//vector<double> seqRead::sub_rate;
gsl_rng* art::gsl_R;
int art::gaussain_mean;
double art::gaussain_sigma;

int main(int argc, char* argv[]){
    wcout << L"==========================================================================="<<endl;
    wcout << L"                           ART (SOLiD version 0.6.2)                       "<<endl; 
    wcout << L" Copyright (c) 2008-2010, Weichun Huang, Jason Myers. All Rights Reserved. "<<endl; 
    wcout << L"==========================================================================="<<endl<<endl;

    bool mask_n=true; 
    short max_num_n=5; 
    int len_ref_id=250;

    //caluate CPUT time
     clock_t start, end;
     double cpu_time_used;
     start = clock();
    
    // Boolean flags showing the state of the run
    bool arg_success = true;
    bool is_pairend_read=false;
    bool mean_flag = false;
    bool sDev_flag = false;
    bool in_flag = false;
    bool out_flag = false;
    bool len_flag = false;
    bool fold_flag = false;
    bool help_flag = false;
    bool show_flag = false;
    bool rate_flag = false;
    bool first_qual = false;
    bool second_qual = false;

    // Command-line Arguments
    string qual_file1 = "";
    string qual_file2 = "";
    char* seq_file = "";
    string out_file_prefix = "";
    int read_len  = 0;
    double x_fold  = 0;
    string num= "";
    int mean = 0;
    double std_dev = 0.0;
    double insRate=0.0001;
    double delRate=0.0001;
    double subRate=0.028;
    double insRate2 = 0.0001;
    double delRate2 = 0.0001;
    double subRate2 = 0.036;

    // Get Command-line arguments
    for(int i=1;i < argc;++i){
        char* arg = argv[i];
            if(!strcmp(arg, "--paired") || !strcmp(arg, "-p") ){
                num="1";
                is_pairend_read=true;
            } else if(!strcmp(arg, "--help") || !strcmp(arg, "-h")){
                arg_success = false;
                help_flag = true;
                break;
            } else if(!strcmp(arg, "--quiet") || !strcmp(arg, "-q")){
                show_flag = true;
            } else if(!strcmp(arg, "--in") || !strcmp(arg, "-i")){
                i++;
                seq_file = argv[i];
                in_flag = true;
            } else if(!strcmp(arg, "--out") || !strcmp(arg, "-o")){
                i++;
                out_file_prefix = argv[i];
                out_flag = true;
            } else if(!strcmp(arg, "--len") || !strcmp(arg, "-l")){
                i++;
                read_len = atoi(argv[i]);
                len_flag = true;
                if(read_len < 0){
                   cerr << "Fatal Error: The read length must be a positive integer." << endl;
                   arg_success = false;
                }
            } else if(!strcmp(arg, "--insRate") || !strcmp(arg, "-ir")){
                i++;
                insRate = atof(argv[i]);
                rate_flag = true;
            } else if(!strcmp(arg, "--delRate") || !strcmp(arg, "-dr")){
                i++;
                delRate = atof(argv[i]);
                rate_flag = true;
            } else if(!strcmp(arg, "--subRate") || !strcmp(arg, "-sr")){
                i++;
                subRate = atof(argv[i]);
                rate_flag = true;
            } else if(!strcmp(arg, "--insRate2") || !strcmp(arg, "-ir2")){
                i++;
                insRate2 = atof(argv[i]);
                rate_flag = true;
            } else if(!strcmp(arg, "--delRate2") || !strcmp(arg, "-dr2")){
                i++;
                delRate2 = atof(argv[i]);
                rate_flag = true;
            } else if(!strcmp(arg, "--subRate2") || !strcmp(arg, "-sr2")){
                i++;
                subRate2 = atof(argv[i]);
                rate_flag = true;
            } else if(!strcmp(arg, "--fcov") || !strcmp(arg, "-f")){
                i++;
                x_fold = atof(argv[i]);
                fold_flag = true;
            } else if(!strcmp(arg, "--mflen") || !strcmp(arg, "-m")){
                i++;
                mean = atoi(argv[i]);
                mean_flag = true;
            } else if(!strcmp(arg, "--sdev") || !strcmp(arg, "-s")){
                i++;
                std_dev = atof(argv[i]);
                sDev_flag = true;
            } else if(!strcmp(arg, "--qprof1") || !strcmp(arg, "-1")){
                i++;
                qual_file1 = argv[i];
		first_qual = true;
            } else if(!strcmp(arg, "--qprof2") || !strcmp(arg, "-2")){
                i++;
                qual_file2 = argv[i];
		second_qual = true;
            } else {
                arg_success = false;
                cerr << "Fatal Error: " << arg << ", is not a valid parameter." << endl;
                break;
            }

    }


    // Make sure the minimum requirements to run were met if the help tag was not given
    if(help_flag){
    } else if(!in_flag || !out_flag || !len_flag || !fold_flag){
        arg_success = false;
        cerr << "Fatal Error: An input-file, output-file prefix, read length, and fold coverage must be"
             << " specified." << endl << endl;
    }

    if(mean_flag && sDev_flag){
        is_pairend_read = true;
        num = "1";
    }

    // Make sure the minimum requirements to run were given for a paired end simulation
    if (is_pairend_read){
        if(mean_flag && sDev_flag){
            art::ini_read_pair_rand(abs(mean),fabs(std_dev));
            if(art::gaussain_mean<=read_len){
                cerr<<"Fatal Error: The read length must be shorter than the mean fragment length specified."
                    <<endl;
                exit(1);
            }
        } else {
            arg_success = false;
            cerr << "Fatal Error: A mean fragment length and a standard deviation must be specified."
                 << endl << endl;
        }
    }

    if(!arg_success){
        wcout << L"------------------------------------USAGE------------------------------------" << endl << endl;
        wcout << L"  -p, --paired       Specify a paired-end read simulation." << endl;
        wcout << L"  -q, --quiet        Turn off end of run summary." << endl;
        wcout << L"  -h, --help         Print out usage information for ART." << endl << endl;
        wcout << L"Parameters should be entered preceeded by the following tags" << endl;
        wcout << L"in any order." << endl << endl;
        wcout << L"  -i, --in           Specify the input file name." << endl;
        wcout << L"  -o, --out          Specify the output file prefix." << endl;
        wcout << L"  -l, --len          Specify the read length to be simulated." << endl;
        wcout << L"  -f, --fcov         Specify the fold coverage to be simulated." << endl;
        wcout << L"  -m, --mflen        Specify the mean fragment length." << endl;
        wcout << L"                     * Only for paired-end simulations." << endl;
        wcout << L"  -s, --sdev         Specify the standard deviation." << endl;
        wcout << L"                     * Only for paired-end simulations." << endl;
        wcout << L"  -ir, --insRate     Specify the first-read insertion rate to be used." << endl;
        wcout << L"                     * The default value is 0.0001." << endl;
        wcout << L"  -dr, --delRate     Specify the first-read deletion rate to be used." << endl;
        wcout << L"                     * The default value is 0.0001." << endl;
        wcout << L"  -sr, --subRate     Specify the first-read substitution rate to be used." << endl;
        wcout << L"                     * The default value is 0.028." << endl;
        wcout << L"  -ir2, --insRate2   Specify the second-read insertion rate to be used." << endl;
        wcout << L"                     * The default value is 0.0001." << endl;
        wcout << L"  -dr2, --delRate2   Specify the second-read deletion rate to be used." << endl;
        wcout << L"                     * The default value is 0.0001." << endl;
        wcout << L"  -sr2, --subRate2   Specify the second-read substitution rate to be used." << endl;
        wcout << L"                     * The default value is 0.036." << endl;
        wcout << L"  -1, --qprof1       Specify the first-read quality profile to be used." << endl;
        wcout << L"                     * If not specified a default profile will be used." << endl;
        wcout << L"  -2, --qprof2       Specify the second-read quality profile to be used." << endl;
        wcout << L"                     * If not specified a default profile will be used." << endl;
        wcout << L"****NOTES****" << endl;
        wcout << L"* To simulate single-end reads one must provide the ART program" << endl;
        wcout << L"  with atleast an input file, output file prefix, the read length," << endl;
        wcout << L"  and the fold coverage." << endl << endl;;
        wcout << L"  Example: art --in reference_DNA.fa --out sim1 --len 35 --fcov 2" << endl << endl;;
        wcout << L"* To simulate paired-end reads one must use the parameters above as" << endl;
        wcout << L"  well as the mean fragment length, standard deviation, and --paired tag." << endl << endl;
        wcout << L"  Example: art --paired --in reference_DNA.fa --out sim2 --len 35 --fcov 2" << endl;
        wcout << L"               --mflen 200 --sdev 3.5" << endl << endl;

        exit(0);
    }
        string seqfasta=out_file_prefix;
        string qualfasta=out_file_prefix;
        string alnfasta=out_file_prefix;
        string fqfile=out_file_prefix;

    if(is_pairend_read){
        seqfasta += "_R3.fa";
        qualfasta += "_R3.qual";
        alnfasta += "_R3.aln";
        fqfile += "_R3.fq";
    } else {
        seqfasta += "_F3.fa";
        qualfasta += "_F3.qual";
        alnfasta += "_F3.aln";
        fqfile += "_F3.fq";
    }

//    ofstream SEQFILE(seqfasta.c_str(),ios::binary);
//    if(!SEQFILE.is_open()) { cout<<"can not open output file: "<<seqfasta<<endl; exit(0); }

//    ofstream QUALFILE(qualfasta.c_str(),ios::binary);
//    if(!QUALFILE.is_open()) { cout<<"can not open output file: "<<qualfasta<<endl; exit(0); }

    ofstream FQFILE(fqfile.c_str(),ios::binary);
    if(!FQFILE.is_open()) { cerr<<"Fatal Error: Can not open output file: "<<fqfile<<endl; exit(0); }

    ofstream ALNFILE(alnfasta.c_str(),ios::binary);
    if(!ALNFILE.is_open()) { cerr<<"Fatal Error: Can not open output file: "<<alnfasta<<endl; exit(0); }

    empdist::set_err_prob();
    empdist qdist;


    vector<short> qual;
    readSeqFile seq_reader(seq_file);
    string id;
    art a_art; 
    seqRead a_read(qual_file1, true);
    a_read.ini_ran_qual();

    if(read_len>a_read.error_profile.size()){
        cerr<<"Error: The read length "<<read_len<<" exceeds the maximum read length:"<<a_read.error_profile.size()<<endl;
        exit(1);
    }
/*
    for(size_t i=0; i<a_read.cal_err_rate_1st.size(); i++){
      a_read.cal_err_rate_1st[i]*=error_scale_factor;
    }
    for(size_t i=0; i<a_read.cal_err_rate_2nd.size(); i++){
      a_read.cal_err_rate_2nd[i]*=error_scale_factor;
    }
*/
    a_read.set_rate(read_len,insRate,2,a_read.ins_rate);
    a_read.set_rate(read_len,delRate,2,a_read.del_rate);
    a_read.set_rate(read_len,subRate,2,a_read.sub_rate);
    string aln_read,aln_ref;
    ostringstream osID;
    int num_seq=0;
    string read_id;

    unsigned long cc_num_read=1;

	string seqfasta2="";
	string qualfasta2="";
	string alnfasta2="";
	string fqfile2="";
    if(is_pairend_read){
        seqfasta2=out_file_prefix+"_F3.fa";
        qualfasta2=out_file_prefix+"_F3.qual";
        alnfasta2=out_file_prefix+"_F3.aln";
        fqfile2=out_file_prefix+"_F3.fq";
//        ofstream SEQFILE2(seqfasta2.c_str(),ios::binary);
//        if(!SEQFILE2.is_open()) { cout<<"can not open output file: "<<seqfasta2<<endl; exit(0); }

//        ofstream QUALFILE2(qualfasta2.c_str(),ios::binary);
//        if(!QUALFILE2.is_open()) { cout<<"can not open output file: "<<qualfasta2<<endl; exit(0); }

        ofstream FQFILE2(fqfile2.c_str(),ios::binary);
        if(!FQFILE2.is_open()) { cerr<<"Fatal Error: Can not open output file: "<<fqfile2<<endl; exit(0); }

        ofstream ALNFILE2(alnfasta2.c_str(),ios::binary);
        if(!ALNFILE2.is_open()) { cerr<<"Fatal Error: Can not open output file: "<<alnfasta2<<endl; exit(0); }
        seqRead a_read_2(qual_file2, false); 
        a_read_2.ini_ran_qual();
/*
        for(size_t i=0; i<a_read_2.cal_err_rate_1st.size(); i++){
          a_read_2.cal_err_rate_1st[i]*=error_scale_factor;
        }
        for(size_t i=0; i<a_read_2.cal_err_rate_2nd.size(); i++){
          a_read_2.cal_err_rate_2nd[i]*=error_scale_factor;
        }
*/
        a_read_2.set_rate(read_len,insRate2,2,a_read.ins_rate);
        a_read_2.set_rate(read_len,delRate2,2,a_read.del_rate);
        a_read_2.set_rate(read_len,subRate2,2,a_read.sub_rate);
        vector<short> qual_2;
        string read_id_2;
        string aln_read_2,aln_ref_2;
        while(seq_reader.next_seq(id,a_art.ref_seq)){ 
//            size_t p1=id.find_first_of(' '); if(p1==string::npos) p1=10; size_t p2=id.find_first_of('\t'); if(p2==string::npos) p2=10;            p1=p1<p2?p1:p2; id=id.substr(0,p1); 
            istringstream isID; isID.str(id); isID>>id; id=id.substr(0,len_ref_id); 
            num_seq++;
            a_art.ini_set(read_len);
            if(mask_n){ 
              a_art.mask_n_region(max_num_n);
            }
            long t_num_read=static_cast<unsigned long>(a_art.ref_seq.size()/read_len*x_fold);
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

                a_read.clear();
                a_read_2.clear();
                //a_art.next_read(a_read);
                //a_art.next_pair_read_indel(a_read, a_read_2); 
                a_art.next_pair_read(a_read, a_read_2); 
                if(mask_n){ 
                  if(a_read.is_plus_strand){
                    if(a_art.masked_pos.count(a_read.bpos)>0 || a_art.masked_pos.count(a_read_2.bpos)>0){
                      t_num_read-=2;
                      continue;
                    }
                  }
                  else{
                    size_t bpos1=a_art.ref_seq.size()-a_read.bpos-read_len;
                    size_t bpos2=a_art.ref_seq.size()-a_read_2.bpos-read_len;
                    if(a_art.masked_pos.count(bpos1)>0 || a_art.masked_pos.count(bpos2)>0){
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
                read_id_2=read_id+"_F3";
                read_id+="_R3";
//print first read
//                SEQFILE<<">"<<read_id<<endl<<a_read.seq_read<<endl; //<<a_read.seq_ref<<endl;
//                QUALFILE<<">"<<read_id<<endl;
//                copy(qual.begin(),qual.end(), ostream_iterator<short>(QUALFILE,"\t"));
//                QUALFILE<<endl;

                FQFILE<<"@"<<read_id<<endl<<'G'<<cs_seq_1st<<endl<<"+"<<endl;
		FQFILE<<"!";
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
		FQFILE2<<"!";
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
            }
        }
//        SEQFILE2.close();
//        QUALFILE2.close();
        FQFILE2.close();
        ALNFILE2.close();
    }
    else{
        while(seq_reader.next_seq(id,a_art.ref_seq)){
            istringstream isID; isID.str(id); isID>>id; id=id.substr(0,len_ref_id); 
            num_seq++;
            a_art.ini_set(read_len);
            long t_num_read=static_cast<unsigned long>(a_art.ref_seq.size()/read_len*x_fold);
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
                a_art.next_read(a_read);
                if(mask_n){ 
                  if(a_read.is_plus_strand){
                    if(a_art.masked_pos.count(a_read.bpos)>0){
                      t_num_read-=1;
                      continue;
                    }
                  }
                  else{
                    size_t bpos=a_art.ref_seq.size()-a_read.bpos-read_len;
                    if(a_art.masked_pos.count(bpos)>0){
                      t_num_read-=1;
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
		FQFILE<<"!";
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
            }
        }
    }

//    SEQFILE.close();
//    QUALFILE.close();
    FQFILE.close();
    ALNFILE.close();

    if(!is_pairend_read){
        wcout << L"                        Single-end Simulation Complete" << endl << endl;
    } else {
        wcout << L"                        Paired-end Simulation Complete" << endl << endl;
    }


    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    wcout<< L"Total CPU time used: "<<cpu_time_used<<endl << endl;

    if(!show_flag){

        wcout << L"Parameters used during run:" << endl << endl;
        wcout << L"\tRead Length:              " << read_len << endl;
        wcout << L"\tFold Coverage:            " << x_fold << "X" << endl;
        if(is_pairend_read){
            wcout << L"\tMean Fragment Length:     " << mean << endl;
            wcout << L"\tStandard Deviation:       " << std_dev << endl;
        }
        if(rate_flag){
            wcout << L"\tFirst Insertion Rate:     " << insRate << endl;
            wcout << L"\tSecond Insertion Rate:    " << insRate2 << endl;
            wcout << L"\tFirst Deletion Rate:      " << delRate << endl;
            wcout << L"\tSecond Deletion Rate:     " << delRate2 << endl;
            wcout << L"\tFirst Substitution Rate:  " << subRate << endl;
            wcout << L"\tSecond Substitution Rate: " << subRate2 << endl;
        }

        wcout << endl << L"Quality Profile(s) Used:" << endl << endl;
        if(is_pairend_read){
            if(first_qual){
                wcout << L"\tFirst Read:  " << qual_file1.c_str()<< endl;
            } else {
                wcout << L"\tFirst Read:  " << "DEFAULT" << endl;
            }
            if(second_qual){
                wcout << L"\tSecond Read: " << qual_file2.c_str()<< endl<< endl;
            } else {
                wcout << L"\tSecond Read: " << "DEFAULT" << endl << endl;
            }
        } else {
            if(first_qual){
                wcout << L"\t" << qual_file1.c_str()<< endl << endl;
            } else {
                wcout << L"\tDEFAULT" << endl << endl;
            }
        }

        wcout << L"Output files created:" << endl << endl;

        if(is_pairend_read){
            wcout << L"  Sequence Files:" << endl;
            wcout << L"\t" << fqfile.c_str() << endl;
            wcout << L"\t" << fqfile2.c_str() << endl << endl;
            wcout << L"  Alignment Files:" << endl;
            wcout << L"\t" << alnfasta.c_str() << endl;
            wcout << L"\t" << alnfasta2.c_str() << endl << endl;;
        } else {
            wcout << L"  Sequence File:" << endl;
            wcout << L"\t" << fqfile.c_str() << endl << endl;
            wcout << L"  Alignment File:" << endl;
            wcout << L"\t" << alnfasta.c_str() << endl << endl;
        }
    }

   return 0;

}

bool art::next_read_indel(seqRead& a_read){
    long pos=(long) floor(r_prob()*valid_region); //pos in [0 ..len-1]   
    int slen =a_read.get_indel(read_len);            
    a_read.is_plus_strand=true;
    if(r_prob()>0.5){
        a_read.is_plus_strand=false;
    }
    if(a_read.is_plus_strand){
        a_read.seq_ref=ref_seq.substr(pos, read_len-slen);
    }
    else{
        a_read.seq_ref=ref_seq_cmp.substr(pos, read_len-slen);
    }
    a_read.bpos=pos;
    a_read.ref2read();
    return true;
}

bool art::next_read(seqRead& a_read){
    long pos=(long) floor(r_prob()*valid_region); //pos in [0 ..len-1]           
    //string seq_ref;
    a_read.is_plus_strand=true;
    if(r_prob()>0.5){
        a_read.is_plus_strand=false;
    }
    if(a_read.is_plus_strand){
        a_read.seq_ref=ref_seq.substr(pos, read_len);
    }
    else{
        a_read.seq_ref=ref_seq_cmp.substr(pos, read_len);
    }
    a_read.bpos=pos;
    a_read.seq_read=a_read.seq_ref;
    return true;
}

bool art::next_pair_read(seqRead& read_1, seqRead& read_2){
    int fragment_len=gaussain_mean+ (int)floor(gsl_ran_gaussian(gsl_R, gaussain_sigma));
    while (fragment_len<read_len || fragment_len>ref_seq.length()){
        fragment_len=gaussain_mean+ (int)floor(gsl_ran_gaussian(gsl_R, gaussain_sigma));
    }
    long pos_1=(long) floor((ref_seq.length()-fragment_len)*r_prob());
    long pos_2=pos_1+fragment_len-read_len;
    bool is_plus_strand=true;
    if(r_prob()>0.5){
        is_plus_strand=false;
    }
    if(is_plus_strand){
        read_1.seq_ref=ref_seq.substr(pos_1, read_len);
        read_2.seq_ref=ref_seq.substr(pos_2, read_len);
    }
    else{
        read_1.seq_ref=ref_seq_cmp.substr(pos_1, read_len);
        read_2.seq_ref=ref_seq_cmp.substr(pos_2, read_len);
    }
    read_1.is_plus_strand=is_plus_strand;
    read_1.bpos=pos_1;
    read_2.is_plus_strand=is_plus_strand;
    read_2.bpos=pos_2;
    return true;
}

bool art::next_pair_read_indel(seqRead& read_1, seqRead& read_2){
    int fragment_len=gaussain_mean+ (int)floor(gsl_ran_gaussian(gsl_R, gaussain_sigma));
    while (fragment_len<read_len || fragment_len>ref_seq.length()){
        fragment_len=gaussain_mean+ (int)floor(gsl_ran_gaussian(gsl_R, gaussain_sigma));
    }
    long pos_1=(long) floor((ref_seq.length()-fragment_len)*r_prob());
    long pos_2=pos_1+fragment_len-read_len;
    int slen_1 =read_1.get_indel(read_len);
    int slen_2 =read_2.get_indel(read_len);   
    bool is_plus_strand=true;
    if(r_prob()>0.5){
        is_plus_strand=false;
    }
    if(is_plus_strand){
        read_1.seq_ref=ref_seq.substr(pos_1, read_len-slen_1);
        read_2.seq_ref=ref_seq.substr(pos_2, read_len-slen_2);
    }
    else{
        read_1.seq_ref=ref_seq_cmp.substr(pos_1, read_len-slen_1);
        read_2.seq_ref=ref_seq_cmp.substr(pos_2, read_len-slen_2);
    }
    read_1.is_plus_strand=is_plus_strand;
    read_1.bpos=pos_1;
    read_1.ref2read();
    read_2.is_plus_strand=is_plus_strand;
    read_2.bpos=pos_2;
    read_2.ref2read();
    //cout<<pos_1<<" a: "<<read_1.seq_read<<endl<<pos_2<<" b: "<<read_2.seq_read<<endl;
    return true;
}
//second read is reverse complemenaty strand 
bool art::next_pair_read_indel_cmp(seqRead& read_1, seqRead& read_2){
    int fragment_len=gaussain_mean+ (int)floor(gsl_ran_gaussian(gsl_R, gaussain_sigma));
    while (fragment_len<read_len || fragment_len>ref_seq.length()){
        fragment_len=gaussain_mean+ (int)floor(gsl_ran_gaussian(gsl_R, gaussain_sigma));
    }
    long pos_1=(long) floor((ref_seq.length()-fragment_len)*r_prob());
    //long pos_2=pos_1+fragment_len-read_len;
    long pos_2=ref_seq.length()-pos_1-fragment_len;
    int slen_1 =read_1.get_indel(read_len);
    int slen_2 =read_2.get_indel(read_len);   
    bool is_plus_strand=true;
    if(r_prob()>0.5){
        is_plus_strand=false;
    }
    if(is_plus_strand){
        read_1.is_plus_strand=true;
        read_1.seq_ref=ref_seq.substr(pos_1, read_len-slen_1); 
        read_2.is_plus_strand=false;
//      pos_2=ref_seq.length()-pos_2-read_len;
        read_2.seq_ref=ref_seq_cmp.substr(pos_2, read_len-slen_2);
    }
    else{
        read_1.is_plus_strand=false;
        read_1.seq_ref=ref_seq_cmp.substr(pos_1, read_len-slen_1);
//      pos_2=ref_seq.length()-pos_2-read_len;
        read_2.is_plus_strand=true;
        read_2.seq_ref=ref_seq.substr(pos_2, read_len-slen_2);
    }
    read_1.bpos=pos_1;
    read_1.ref2read();
    read_2.bpos=pos_2;
    read_2.ref2read();
    //cout<<pos_1<<" a: "<<read_1.seq_read<<endl<<pos_2<<" b: "<<read_2.seq_read<<endl;
    return true;
}
//

