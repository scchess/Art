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

#include <cmath>
#include <cstdlib>
#include <ctime>
#include "read_profile.h"

void read_profile::default_profile(bool is_titanium){
    istringstream  distss;
    if(is_titanium){
	    distss.str(t_QUAL_DIST_1BASE);
	    read_emp_dist(distss,qual_1stbase_dist);
	    distss.clear();
	    distss.str(t_QUAL_DIST_MC);
	    read_emp_dist(distss,qual_dist_mc, 3);
	    distss.clear();
	    distss.str(t_READ_LEN_DIST);
	    read_emp_dist(distss,read_len_dist);
    }
    else{
	    distss.str(QUAL_DIST_1BASE);
	    read_emp_dist(distss,qual_1stbase_dist);
	    distss.clear();
	    distss.str(QUAL_DIST_MC);
	    read_emp_dist(distss,qual_dist_mc, 3);
	    distss.clear();
	    distss.str(READ_LEN_DIST);
	    read_emp_dist(distss,read_len_dist);
    }
    const int NUM=10;
    double Undercall[10]={0.0000000000, 0.001096668, 0.004562394, 0.01875503, 0.05466385, 0.14088050, 0.26817640, 0.44392523, 0.54782609, 0.69696970};
    double Overcall[]={0.0005165388, 0.000858758, 0.003219707, 0.01025744, 0.01770155, 0.01677149, 0.02264601, 0.01401869, 0.02608696, 0.03030303};
    double OverErroRate[]={0.0005165388, 0.001955426, 0.007782101, 0.02901247, 0.07236540, 0.15765199, 0.29082241, 0.45794393, 0.57391304, 0.72727273};
    under_call_error.assign(Undercall,Undercall+NUM);
    over_call_error.assign(Overcall,Overcall+NUM);
    total_error.assign(OverErroRate,OverErroRate+NUM);
    cyc_base[0]='T'; cyc_base[1]='A'; cyc_base[2]='C'; cyc_base[3]='G';
}

void read_profile::user_profile(string qual_profile, string mc_qual_profile, string indel_erro_profile, string length_profile){
    read_emp_dist(qual_profile,qual_1stbase_dist);
    read_emp_dist(mc_qual_profile,qual_dist_mc, 3);
    read_indel_error(indel_erro_profile);
    read_emp_dist(length_profile,read_len_dist);
    size_t k=under_call_error.size();
    if(k>over_call_error.size()) k=over_call_error.size();
    total_error.resize(k);
    for(size_t i=0; i<k; i++){
	    total_error[i]=over_call_error[i]+under_call_error[i];
    }
    cyc_base[0]='T'; cyc_base[1]='A'; cyc_base[2]='C'; cyc_base[3]='G';
}

void read_profile::read_indel_error(string error_profile_file){
	//if file is blank, use the default error profile
    ifstream input(error_profile_file.c_str());
    if(!input){
		const int NUM=10;
	       	double Undercall[10]={0.0000000000, 0.001096668, 0.004562394, 0.01875503, 0.05466385, 0.14088050, 0.26817640, 0.44392523, 0.54782609, 0.69696970};
	       	double Overcall[]={0.0005165388, 0.000858758, 0.003219707, 0.01025744, 0.01770155, 0.01677149, 0.02264601, 0.01401869, 0.02608696, 0.03030303};
	       	double OverErroRate[]={0.0005165388, 0.001955426, 0.007782101, 0.02901247, 0.07236540, 0.15765199, 0.29082241, 0.45794393, 0.57391304, 0.72727273};
	       	under_call_error.assign(Undercall,Undercall+NUM);
	       	over_call_error.assign(Overcall,Overcall+NUM);
	       	total_error.assign(OverErroRate,OverErroRate+NUM);
//	       	cerr<<"Note: "<<error_profile_file<<" does not exist. The default built-in indel_error_profile is used for the simulation"<<endl;
		return;
    }
    bool is_UC_profile=false;
    bool is_OC_profile=false;
    while(!input.eof()){
	string profile;
        string aLine;
        getline(input, aLine);
        if(aLine.length()==0) continue;
        if(aLine[0]=='#') continue;
        if(aLine[0]=='*') break;
        istringstream ss(aLine);
        ss>>profile;
        if(profile=="UNDERCALL"){
	       	double t_double;
	       	while (ss >> t_double){ under_call_error.push_back(t_double); }
		if(under_call_error.size()>0){ is_UC_profile=true; }
        }
	else if(profile =="OVERCALL"){
	       	double t_double;
	       	while (ss >> t_double){ over_call_error.push_back(t_double); }
		if(over_call_error.size()>0){ is_OC_profile=true;}
	}
	
    }
    if(!is_UC_profile){
	    cerr<<"Error: there is no the under-call error profile in the input file:"<<error_profile_file<<endl;
	    exit(1);
    }
    if(!is_OC_profile){
	    cerr<<"Error: there is no the over-call error profile in the input file:"<<error_profile_file<<endl;
	    exit(1);
    }

}

bool read_profile::get_read_qual_fast(vector<int>& lenVec, unsigned long b_pos, string& read_seq, string& seq, vector<string>&aln, vector<short>& qual, short& cyc, short cycle){
    if(cyc>=cycle){
	    //cerr<<"initial flow cycle number is larger than maximum cycle number"<<endl;
	    return false;
	    //exit(1);
    }
    if(read_seq.length()==0){
	    //cerr<<"len of read ==0"<<endl;
	    return false;
    }
    if(qual_1stbase_dist.size()==0) return false;
    aln.clear();
    aln.resize(2);
    aln[0].reserve(read_seq.size()+10);
    aln[1].reserve(read_seq.size()+10);
    seq.clear();
    seq.reserve(read_seq.size()+10);
    qual.reserve(read_seq.size()+10);
    size_t len=read_seq.length();
    char aBase, orgBase, cycB;
    short polymer_len;
    for(size_t i=0, k=0; i<len; k++){
	    if(k>3){ k=0; cyc++; } 
	    //stop if max_cycle reached
	    if(cyc>cycle){
		    read_seq=read_seq.substr(0,i);
		    break;
	    }
	    cycB=cyc_base[k];
	    aBase=read_seq[i];
	    if(cycB!=aBase){ continue; }

	    orgBase=aBase;
        polymer_len=lenVec[b_pos+i];
        i+=polymer_len;
       	//last homolymer (trucated)
	if(i>len){
	    int kk=i-len; i=len; 
	    polymer_len-=kk;
       	} 
        //polymer length >that defined in distribution
        short qvalue=60;
        if(polymer_len>qual_1stbase_dist.size()){
            short min_qvalue=30;
            int org_len=polymer_len; 
            while(polymer_len>qual_1stbase_dist.size()){
                qual.push_back(qvalue);
                qvalue-=2;
                if(qvalue<=min_qvalue) qvalue=min_qvalue;
                polymer_len--;
            }
            aln[0].append(org_len-polymer_len, orgBase);
            aln[1].append(org_len-polymer_len, aBase);
            seq.append(org_len-polymer_len,aBase);
        }

        //assign qual of first base
        unsigned int cumCC;
        map <unsigned int, unsigned short>::iterator it;
        cumCC=(unsigned int) ceil(r_prob()*max_dist_number);
        it=qual_1stbase_dist[polymer_len-1].lower_bound(cumCC);
        //qual.push_back(it->second);
        qvalue=it->second;

        //add substitution error for single base only
        if(polymer_len==1){
            if(r_prob()<prob_err[qvalue]){
                char achar=aBase;
                while(aBase==achar){ achar=rand_base(); }
                aBase=achar;
            }
        }

        //add indel error
        short ins_len=0, del_len=0; //
        double p_ran=r_prob();
        if(polymer_len>total_error.size()){ //deal with undefined in error distribution 
            double p_del=total_error[total_error.size()-1];
            p_del=p_del+(1-p_del)*polymer_len/(polymer_len+100); //p deletion p_del->1 as ploymer len ->larger
            double p_ins=(1-p_del)/10; //p insetion  
            if(p_ran<(p_del+p_ins)){
                if(p_ran<p_ins){
                    do{
                        ins_len+=1;
                        p_ins=p_ins*exp((double) -2*ins_len); //exponetial decay 
                    } while(polymer_len>ins_len && r_prob()<p_ins);
                }
                else{
                    do{
                        del_len+=1;
                        p_del=p_del*exp((double) -2*del_len); //exponetial decay 
                    } while(polymer_len>del_len && r_prob()<p_del);
                }
            }
        }
        else if(p_ran<total_error[polymer_len-1]){
            double p_ins=under_call_error[polymer_len-1];
            double p_del=over_call_error[polymer_len-1];
            if(p_ran<p_ins){//insertion
                do{
                    ins_len+=1;
                    p_ins=p_ins*exp((double) -2*ins_len); //exponetial decay 
                } while(polymer_len>ins_len && r_prob()<p_ins);
            }
            else{ //deletion
                do{
                    del_len+=1;
                    p_del=p_del*exp((double) -2*del_len); //exponetial decay 
                } while(polymer_len>del_len && r_prob()<p_del);
            }
        }

        if(ins_len>0){
            aln[0].append(polymer_len, orgBase);
            aln[0].append(ins_len, '-');
            aln[1].append(polymer_len+ins_len, aBase);
            seq.append(polymer_len+ins_len,aBase);
        }
        else if(del_len>0){
            aln[0].append(polymer_len, orgBase);
            aln[1].append(polymer_len-del_len, aBase);
            aln[1].append(del_len, '-');
            seq.append(polymer_len-del_len,aBase);
        }
        else{
            aln[0].append(polymer_len, orgBase);
            aln[1].append(polymer_len, aBase);
            seq.append(polymer_len, aBase);
        }

        //the polymer_len after indel
        polymer_len+=ins_len-del_len;
        //check whether it is delete completely, otherwise add base quality 
        if(polymer_len>0) qual.push_back(qvalue);
        for(short k=1; k<polymer_len; k++){
            ostringstream osID;
            osID<<polymer_len<<"-"<<k<<"-"<<qvalue;
            string idd=osID.str();
            cumCC=(unsigned int) ceil(r_prob()*max_dist_number);
            if(qual_dist_mc.count(idd)>0){
                it=qual_dist_mc[idd].lower_bound(cumCC);
                //if(it!=qual_dist[i].end()){
                qvalue=it->second;
            }
            else qvalue=1;
            qual.push_back(qvalue);
        }

    }
       	return true;
}



bool read_profile::get_read_qual(string& read_seq, string& seq, vector<string>&aln, vector<short>& qual, short cycle){
    if(read_seq.length()==0) return false;
    if(qual_1stbase_dist.size()==0) return false;
    aln.clear();
    aln.resize(2);
    aln[0].reserve(read_seq.size()+10);
    aln[1].reserve(read_seq.size()+10);
    seq.clear();
    seq.reserve(read_seq.size()+10);
    qual.reserve(read_seq.size()+10);
    short cyc=0;
    size_t len=read_seq.length();
    char aBase, orgBase, cycB;
    short polymer_len;
    for(size_t i=0, k=0; i<len; k++){
	    if(k>3){ k=0; cyc++; } 
	    //stop if max_cycle reached
	    if(cyc>cycle){
		    read_seq=read_seq.substr(0,i);
		    break;
	    }
	    cycB=cyc_base[k];
	    aBase=read_seq[i];
	    if(cycB!=aBase){ continue; }
        orgBase=aBase;
        polymer_len=0;
        while(i<len && aBase==read_seq[i]){
            i++;
            polymer_len++;
        }

        //polymer length >that defined in distribution
        short qvalue=60;
        if(polymer_len>qual_1stbase_dist.size()){
            short min_qvalue=30;
            int org_len=polymer_len; 
            while(polymer_len>qual_1stbase_dist.size()){
                qual.push_back(qvalue);
                qvalue-=2;
                if(qvalue<=min_qvalue) qvalue=min_qvalue;
                polymer_len--;
            }
            aln[0].append(org_len-polymer_len, orgBase);
            aln[1].append(org_len-polymer_len, aBase);
        }

        //assign qual of first base
        unsigned int cumCC;
        map <unsigned int, unsigned short>::iterator it;
        cumCC=(unsigned int) ceil(r_prob()*max_dist_number);
        it=qual_1stbase_dist[polymer_len-1].lower_bound(cumCC);
        //qual.push_back(it->second);
        qvalue=it->second;

        //add substitution error for signal base only
        if(polymer_len==1){
            if(r_prob()<prob_err[qvalue]){
                char achar=aBase;
                while(aBase==achar){ achar=rand_base(); }
                aBase=achar;
            }
        }

        //add indel error
        short ins_len=0, del_len=0; //
        double p_ran=r_prob();
        if(polymer_len>total_error.size()){ //deal with undefined in error distribution 
            double p_del=total_error[total_error.size()-1];
            p_del=p_del+(1-p_del)*polymer_len/(polymer_len+100); //p deletion p_del->1 as ploymer len ->larger
            double p_ins=(1-p_del)/10; //p insetion  
            if(p_ran<(p_del+p_ins)){
                if(p_ran<p_ins){
                    do{
                        ins_len+=1;
                        p_ins=p_ins*exp((double) -2*ins_len); //exponetial decay 
                    } while(polymer_len>ins_len && r_prob()<p_ins);
                }
                else{
                    do{
                        del_len+=1;
                        p_del=p_del*exp((double) -2*del_len); //exponetial decay 
                    } while(polymer_len>del_len && r_prob()<p_del);
                }
            }
        }
        else if(p_ran<total_error[polymer_len-1]){
            double p_ins=under_call_error[polymer_len-1];
            double p_del=over_call_error[polymer_len-1];
            if(p_ran<p_ins){//insertion
                do{
                    ins_len+=1;
                    p_ins=p_ins*exp((double) -2*ins_len); //exponetial decay 
                } while(polymer_len>ins_len && r_prob()<p_ins);
            }
            else{ //deletion
                do{
                    del_len+=1;
                    p_del=p_del*exp((double) -2*del_len); //exponetial decay 
                } while(polymer_len>del_len && r_prob()<p_del);
            }
        }

        if(ins_len>0){
            aln[0].append(polymer_len, orgBase);
            aln[0].append(ins_len, '-');
            aln[1].append(polymer_len+ins_len, aBase);
            seq.append(polymer_len+ins_len,aBase);
        }
        else if(del_len>0){
            aln[0].append(polymer_len, orgBase);
            aln[1].append(polymer_len-del_len, aBase);
            aln[1].append(del_len, '-');
            seq.append(polymer_len-del_len,aBase);
        }
        else{
            aln[0].append(polymer_len, orgBase);
            aln[1].append(polymer_len, aBase);
            seq.append(polymer_len, aBase);
        }

        //the polymer_len after indel
        polymer_len+=ins_len-del_len;
        //check whether it is delete completely, otherwise add base quality 
        if(polymer_len>0) qual.push_back(qvalue);
        for(short k=1; k<polymer_len; k++){
            ostringstream osID;
            osID<<polymer_len<<"-"<<k<<"-"<<qvalue;
            string idd=osID.str();
            cumCC=(unsigned int) ceil(r_prob()*max_dist_number);
            if(qual_dist_mc.count(idd)>0){
                it=qual_dist_mc[idd].lower_bound(cumCC);
                //if(it!=qual_dist[i].end()){
                qvalue=it->second;
            }
            else qvalue=1;
            qual.push_back(qvalue);
        }

    }
       	return true;
}


bool read_profile::get_read_qual(vector< map <unsigned int, unsigned short> >&qual_dist, vector<short>& read_qual, int len){
    if((int)qual_dist.size()<len) return false;
    unsigned int cumCC;
    map <unsigned int, unsigned short>::iterator it;
    for(int i=0; i<len; i++){
        cumCC=(unsigned int) ceil(r_prob()*max_dist_number);
        it=qual_dist[i].lower_bound(cumCC);
        //if(it!=qual_dist[i].end()){
        read_qual.push_back(it->second);
        //}
    }
    return true;
}

read_profile::~read_profile(){
    qual_1stbase_dist.clear();
    qual_dist_mc.clear();
}

bool read_profile::read_emp_dist(string infile, vector< map <unsigned int, unsigned short> >& qual_dist){
    ifstream distss(infile.c_str());
    if(!distss){
        cerr<<"Error: cannot open the profile file "<<infile<<endl;
        exit(1);
    }
    return read_emp_dist(distss, qual_dist);
}


bool read_profile::read_emp_dist(istream& input, vector< map <unsigned int, unsigned short> >& qual_dist){
    int linenum=1; //one-based index
    int read_pos;
    while(!input.eof()){
        string aLine;
        getline(input, aLine);
        if(aLine.length()==0) continue;
        if(aLine[0]=='#') continue;
        if(aLine[0]=='*') break;
        istringstream ss(aLine);
        ss>>read_pos;
        if(read_pos!=linenum){
            cerr<<read_pos<<"\t"<<linenum<<endl;
            cerr<<"Fatal error (1): wrong format of input distribution at: "<<endl;
            cerr<<aLine<<endl;
            exit(1);
        }
        unsigned short t_int;
        vector<unsigned short> qual;
        while (ss >> t_int){ qual.push_back(t_int); }
        getline(input, aLine);
        ss.clear();ss.str(aLine);
        ss>>read_pos;
        if(read_pos!=linenum){
            cerr<<"Fatal error (2): wrong format of input distribution"<<endl;
            exit(2);
        }
        unsigned long t_uint;
        vector<unsigned long> count;
        if(IS_CUM_DIST){  while (ss >> t_uint){  count.push_back(t_uint); } }
        else {
            int index=0;
            while (ss >> t_uint){
                count.push_back(t_uint);
                if(index>0) count[index]+=count[index-1];
                index++;
            }
        }
        if(count.size()!=qual.size()){
            cerr<<"Fatal error (3): wrong format of input distribution"<<endl;
            exit(3);
        }
        double denom=count[count.size()-1]/(double)max_dist_number;
        map<unsigned int, unsigned short> dist;
        for(size_t i=0; i<count.size(); i++){
            unsigned int cc=(unsigned int)ceil(count[i]/denom);           
            dist[cc]=qual[i];
        }
        if(dist.size()>0){
            linenum++;
            qual_dist.push_back(dist);
        }
    }
    //max_len_read=linenum;
    if(linenum==0) return false;
    return true;
}


bool read_profile::read_emp_dist(string infile, map<string, map <unsigned int, unsigned short> >& qual_dist, short num_id_field){
    ifstream distss(infile.c_str());
    if(!distss){
        cerr<<"Error: cannot open the distribution file "<<infile<<endl;
        exit(1); 
    }
    return read_emp_dist(distss, qual_dist, num_id_field);
}

bool read_profile::read_emp_dist(istream& input, map<string, map <unsigned int, unsigned short> >& qual_dist, short num_id_field){
    int linenum=0;
    if(num_id_field<=0) return false;
    while(!input.eof()){
        string aLine;
        string pos_id;
        getline(input, aLine);
        if(aLine.length()==0) continue;
        if(aLine[0]=='#') continue;
        if(aLine[0]=='*') break;
        istringstream ss(aLine);
        string tmp_id;
        for(short cc=0; cc<num_id_field; cc++){
            ss>>tmp_id; 
            if(cc==0) pos_id=tmp_id;
            else pos_id+="-"+tmp_id;
        }
        //if(pos_id!=linenum){
        //    cerr<<"Fatal error (1): wrong format of input distribution"<<endl;
        //    exit(1);
        //}
        unsigned short t_int;
        vector<unsigned short> qual;
        while (ss >> t_int){ qual.push_back(t_int); }
        getline(input, aLine);
        ss.clear();ss.str(aLine);
        string id2;
        for(short cc=0; cc<num_id_field; cc++){
            ss>>tmp_id;
            if(cc==0) id2=tmp_id;
            else id2+="-"+tmp_id;
        }
        if(id2!=pos_id){
            cerr<<"Fatal error (2): wrong format of input distribution"<<endl;
            exit(2);
        }
        unsigned long t_uint;
        vector<unsigned long> count;
        if(IS_CUM_DIST){  while (ss >> t_uint){  count.push_back(t_uint); } }
        else {
            int index=0;
            while (ss >> t_uint){
                count.push_back(t_uint);
                if(index>0) count[index]+=count[index-1];
                index++;
            }
        }
        if(count.size()!=qual.size()){
            cerr<<"Fatal error (3): wrong format of input distribution"<<endl;
            exit(3);
        }
        double denom=count[count.size()-1]/(double)max_dist_number;
        map<unsigned int, unsigned short> dist;
        for(size_t i=0; i<count.size(); i++){
            unsigned int cc=(unsigned int)ceil(count[i]/denom);           
            dist[cc]=qual[i];
        }
        if(dist.size()>0){
            linenum++;
            qual_dist[pos_id]=dist;
        }
    }
    //max_len_read=linenum;
    if(linenum==0) return false;
    return true;
}

bool read_profile::read_emp_dist(istream& input, map <unsigned int, unsigned int>& dist){
    while(!input.eof()){
        string aLine;
        getline(input, aLine);
        if(aLine.length()==0) continue;
        if(aLine[0]=='#') continue;
        istringstream ss(aLine);
        unsigned int t_int;
        vector<unsigned int> value;
        while (ss >> t_int){ value.push_back(t_int); }
        getline(input, aLine);
        ss.clear();ss.str(aLine);
        vector<unsigned int> count;
        if(IS_CUM_DIST){  while (ss >> t_int){  count.push_back(t_int); } }
        else {
            int index=0;
            while (ss >> t_int){
                count.push_back(t_int);
                if(index>0) count[index]+=count[index-1];
                index++;
            }
        }
        if(count.size()!=value.size()){
            cerr<<"Fatal error (3): wrong format of input distribution"<<endl;
            exit(3);
        }
        double denom=count[count.size()-1]/(double)max_dist_number;
        for(size_t i=0; i<count.size(); i++){
            unsigned int cc=(unsigned int)ceil(count[i]/denom); 
            dist[cc]=value[i];
        }
        break;
    }
    if(dist.size()==0) return false;
    return true;
}

bool read_profile::read_emp_dist(string infile, map <unsigned int, unsigned int>& dist){
    ifstream distss(infile.c_str());
    if(!distss){
        cerr<<"Error: cannot open the distribution file "<<infile<<endl;
        exit(1); 
    }
    return read_emp_dist(distss, dist);
}
