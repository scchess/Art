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

#include <algorithm>
#include "read_SOLiD.h"

SOLiDread::SOLiDread(string profile_file){
       	double error_1st []={};
       	double error_2nd[]={};
       	//color space to sequence
	cs2seq['a']['0'] = 'a';
       	cs2seq['a']['1'] = 'c';
       	cs2seq['a']['2'] = 'g';
       	cs2seq['a']['3'] = 't';
        cs2seq['c']['0'] = 'c';
        cs2seq['c']['1'] = 'a';
        cs2seq['c']['2'] = 't';
        cs2seq['c']['3'] = 'g';
        cs2seq['g']['0'] = 'g';
        cs2seq['g']['1'] = 't';
        cs2seq['g']['2'] = 'a';
        cs2seq['g']['3'] = 'c';
        cs2seq['t']['0'] = 't';
        cs2seq['t']['1'] = 'g';
        cs2seq['t']['2'] = 'c';
        cs2seq['t']['3'] = 'a';
        //sequence to color space
        seq2cs['a']['a'] = '0';
        seq2cs['a']['c'] = '1';
        seq2cs['a']['g'] = '2';
        seq2cs['a']['t'] = '3';
        seq2cs['c']['c'] = '0';
        seq2cs['c']['a'] = '1';
        seq2cs['c']['t'] = '2';
        seq2cs['c']['g'] = '3';
        seq2cs['g']['g'] = '0';
        seq2cs['g']['t'] = '1';
        seq2cs['g']['a'] = '2';
        seq2cs['g']['c'] = '3';
        seq2cs['t']['t'] = '0';
        seq2cs['t']['g'] = '1';
        seq2cs['t']['c'] = '2';
        seq2cs['t']['a'] = '3';
	if(profile_file.empty()){
	       	istringstream  in_profile;
	       	in_profile.str(abSOLiD_PROFILE);
	       	read_error_profile(in_profile);
	       	in_profile.clear();
	}
	else if(profile_file=="pseudo"){
	       	istringstream  in_profile;
	       	in_profile.str(pseudo_PROFILE);
	       	read_error_profile(in_profile);
	       	in_profile.clear();
	}
	else{
	       	read_error_profile(profile_file);
	}
};


int SOLiDread::get_indel(int read_len){
    //if(ins_rate.size()>=read_len) {cerr<<"fatal error\n";  exit(1)};
    int ins_len=0, del_len=0;
    //deletion
//  char base_type[4]={'a','c','g','t'};
    char base_type[4]={'0','1','2','3'}; //color space
    for(int i=(int)del_rate.size()-1; i>=0; i--){
        if(del_rate[i]>=r_prob()){
            del_len=i+1;
            for(int j=i; j>=0;){
                int pos=(int) floor((read_len-1)*r_prob()); //invalid deletion positions: 0 or read_len-1 as these deletions are not shown in reads
                if(pos==0) continue;
                if(indel.count(pos)==0){
                    indel[pos]='-';
                    j--;
                }
            }
            break;
        }
    }

    for(int i=ins_rate.size()-1; i>=0; i--){
        if(ins_rate[i]>=r_prob()){
            ins_len=i+1;
            for(int j=i; j>=0;){
                int pos=(int) floor(r_prob()*read_len);
                if(indel.count(pos)==0){
                    short base=(short)ceil(r_prob()*4);
		    indel[pos]=base_type[base-1];
//                    switch(base){
//                                case 1:
//                                    indel[pos]='a';   break;
//                                case 2:
//                                    indel[pos]='c';   break;
//                                case 3:
//                                    indel[pos]='g';   break;
//                                case 4:
//                                    indel[pos]='t';  
//                    }
                    j--;
                }
            }
            break;
        }
    }
    return (ins_len-del_len);
};

void SOLiDread::ref2read(){
    if(indel.size()==0){
        seq_read=seq_ref;
        return;
    }
    seq_read.clear();
    int k=0;
    for(size_t i=0; i<seq_ref.size();){
        //cout<<i<<"\t"<<k<<endl;
        if(indel.count(k)==0){
            seq_read.push_back(seq_ref[i]); i++; k++; 
        }
        else if(indel[k]=='-'){
            i++;k++;
        }
        else{
            seq_read.push_back(indel[k]); k++;
        }
    }
    while(indel.count(k)>0){
        seq_read.push_back(indel[k]);
        k++;
    }
}

bool SOLiDread::read_error_profile(string infile_profile){
       	ifstream input(infile_profile.c_str(), ifstream::in);
       	if(!input) { cerr<<"can not open error profile file: "<<infile_profile<<endl; exit(1); }
	else{
	       	read_error_profile(input);
	}
       	return false;
};

//read position and color dependent error profile for paired-end reads
bool SOLiDread::read_error_profile(istream& input){
       	//1	A	C	0.001491472	0.003540099
	int last_read_pos=1;
       	int read_pos;
       	map <char, base_quality> aMAP;
       	short k=0;
       	base_quality bq;
       	char cur_base;
       	while(!input.eof()){
	       	string aLine;
	       	getline(input, aLine);
	       	if(aLine.length()==0) continue;
	       	if(aLine[0]=='#') continue;
	       	if(aLine[0]=='*') break;
	       	istringstream ss(aLine);
	       	ss>>read_pos;
	       	if(last_read_pos<read_pos){
		       	error_profile.push_back(aMAP);
		       	aMAP.clear();
		       	last_read_pos=read_pos;
	       	}
//cerr<<read_pos<<"\t"<<error_profile.size()<<endl; 
	       	if(read_pos!=error_profile.size() && (read_pos-1)!=(int)error_profile.size()){
		       	cerr<<"Fatal error (1): wrong format of error profile"<<endl;
		       	cerr<<aLine<<endl;
		       	exit(1);
	       	} 
		char right_b, wrong_b; 
		double error_rate[2]={0, 0}; //prob_error_1st,  prob_error_2nd,  quality_1st,  quality_2nd,  
		ss >> right_b; right_b=tolower(right_b);
	       	ss >> wrong_b; wrong_b=tolower(wrong_b);
	       	ss >> error_rate[0];
	       	if (! (ss >> error_rate[1])){
		       	cerr<<"Fatal error (2): wrong format of error profile"<<endl;
		       	cerr<<aLine<<endl;
		       	exit(1);
	       	}
		if(right_b == wrong_b){
		       	cerr<<"Fatal error (3): wrong format input color base (correct color=wrong color)"<<endl;
		       	cerr<<"\t"<<aLine<<endl;
		       	exit(1);
	       	}
	       	if(k==0){
		       	cur_base=right_b;
	       	}
	       	else{
		       	if(cur_base!=right_b){
			       	cerr<<"Fatal error (3): number of error type is wrong"<<endl;
			       	cerr<<"\t"<<aLine<<endl;
			       	exit(1);
		       	} 
			//cumulative error rate;
			error_rate[0]+=bq.error_prob[0][k-1];
			error_rate[1]+=bq.error_prob[1][k-1];
		}

		bq.error_color[k] = wrong_b;
		bq.error_prob[0][k]=error_rate[0];
		bq.error_prob[1][k]=error_rate[1];

		if(k==2){
		       	aMAP[right_b]=bq;
			aMAP[right_b].qual[0]=-10*log10(bq.error_prob[0][2]);
			//aMAP[right_b][1].prob=bq[1].error_prob[2];
			aMAP[right_b].qual[1]=-10*log10(bq.error_prob[1][2]);
		       	k=0;
	       	}
	       	else{
		       	k++;
	       	}
       	}
       	if(error_profile.size()>0) return true;
       	else return false;
}

void SOLiDread::ini_ran_qual(unsigned int r_seed){
       	gsl_rng_default_seed=r_seed;
       	const gsl_rng_type *rndT=gsl_rng_default;
       	gsl_p_1st.push_back(gsl_rng_alloc(rndT));
       	gsl_p_1st.push_back(gsl_rng_alloc(rndT));
}

void SOLiDread::ini_ran_qual(){
       	gsl_rng_default_seed=(unsigned int)time(NULL);
       	const gsl_rng_type *rndT=gsl_rng_default;
       	gsl_p_1st.push_back(gsl_rng_alloc(rndT));
       	gsl_p_1st.push_back(gsl_rng_alloc(rndT));
} 

void SOLiDread::set_rate(int read_len, double p, int max_num, vector <double>& rate){
       	rate.resize(max_num);
       	for(size_t i=1; i<=max_num; i++){
	       	rate[i-1]= gsl_cdf_binomial_Q(i, p, read_len);
       	}
}

//convert base-space to color space, and incorporate sequencing errors
void SOLiDread::convert_seq2cs(string& c_space, vector<short>& qual, map<int,char>& error_pos, bool is_F3){ 
	int len=seq_ref.length();
	if(len==0) return;
	error_pos.clear();
	short n=0;
	if(!is_F3) n=1;
	c_space.resize(len);
	qual.resize(len);
	//first reads always start with 'T', while 2nd reads start with 'G'
	char lastchar='t';
	if(!is_F3) lastchar='g';
	char achar=tolower(seq_ref[0]); 
	char nchar=seq2cs[lastchar][achar];
       	qual[0]=gsl_ran_poisson(gsl_p_1st[n], error_profile[0][achar].qual[n]);
       	if(qual[0]>max_qual_value) qual[0]=max_qual_value;
       	double p=r_prob();
       	if(p < prob_err[qual[0]]){ 
		p *= prob_err[qual[0]]/error_profile[0][achar].error_prob[n][2]; //scale p to match orignial error probability 
		if(p <error_profile[0][achar].error_prob[n][0]){
		       	c_space[0]=error_profile[0][achar].error_color[0];
	       	}
	       	else if(p < error_profile[0][achar].error_prob[n][1]){
		       	c_space[0]=error_profile[0][achar].error_color[1];
	       	}
	       	else{
		       	c_space[0]=error_profile[0][achar].error_color[2];
	       	}
	       	c_space[0]=seq2cs[lastchar][c_space[0]];
	       	error_pos[0]=nchar;
       	}
       	else{
	       	c_space[0]=nchar;
       	}

	for(int i=1; i<len; i++){ 
		achar=tolower(seq_ref[i]); 
		//if(i>0) lastchar=tolower(seq_ref[i-1]);
		lastchar=tolower(seq_ref[i-1]);
		nchar=seq2cs[lastchar][achar];
		qual[i]=gsl_ran_poisson(gsl_p_1st[n], error_profile[i][nchar].qual[n]);
		if(qual[i]>max_qual_value) qual[i]=max_qual_value;
		double p=r_prob();
	       	if(p < prob_err[qual[i]]){ 
			p *= prob_err[qual[i]]/error_profile[i][nchar].error_prob[n][2]; //scale p to match orignial error probability 
			if(p <error_profile[i][nchar].error_prob[n][0]){
			      	c_space[i]=error_profile[i][nchar].error_color[0];
			}
			else if(p < error_profile[i][nchar].error_prob[n][1]){
			       	c_space[i]=error_profile[i][nchar].error_color[1];
			}
			else{
			       	c_space[i]=error_profile[i][nchar].error_color[2];
			}
		       	error_pos[i]=nchar;
		}
		else{
		       	c_space[i]=nchar;
		}
	}
}

//convert base-space to color space
void SOLiDread::convert_seq2cs(string& b_space, string& c_space){ 
	int len=b_space.length();
	if(len==0) return;
	c_space.resize(len);
	c_space[0]=seq2cs['t'][tolower(b_space[0])];
       	for(int i=1; i<len; i++){
		c_space[i]=seq2cs[tolower(b_space[i-1])][tolower(b_space[i])];
       	}
}

//convert color space to base space
void SOLiDread::convert_cs2seq(string& c_space,  string& b_space, char base1st){ 
	int len=c_space.length();
	if(len==0) return;
	b_space.resize(len);
	b_space[0]=cs2seq[base1st][c_space[0]];
       	for(int i=1; i<len; i++){
		b_space[i]=cs2seq[b_space[i-1]][c_space[i]];
       	}
}

char SOLiDread::rand_base(){
       	short base=(short)ceil(r_prob()*4);
        switch(base){
            case 1:
                return 'a';
            case 2:
                return 'c';
            case 3:
                return 'g';
            case 4:
                return 't';  
            default:
		cerr<<"system error.";
                exit(1); 
        }
}

bool SOLiDread::get_aln(string& aln_read, string& aln_ref){
       	if(indel.size()==0) return false;
        map<int,char,less<int> >::iterator it;
        aln_read=seq_read; aln_ref=seq_ref;
        for(it=indel.begin(); it!=indel.end(); it++){
            if(it->second!='-'){
                aln_ref.insert(it->first,1,'-');
            }
            else{
                aln_read.insert(it->first,1,'-');
            }
        }
        return true;
}

//based on based on calibrated position-depended error rates
int SOLiDread::add_calib_error_1st(){
       	int num=0;
        for(size_t i=0; i<seq_read.size(); i++){
            if(r_prob()<cal_err_rate_1st[i]){
                char achar=seq_read[i];
                while(seq_read[i]==achar){ achar=rand_base(); }
                seq_read[i]=achar;
                substitution[i]=achar;
                num++;
            }
        }
        return num;
}

int SOLiDread::add_calib_error_2nd(){
        int num=0;
        for(size_t i=0; i<seq_read.size(); i++){
            if(r_prob()<cal_err_rate_2nd[i]){
                char achar=seq_read[i];
                while(seq_read[i]==achar){ achar=rand_base(); }
                seq_read[i]=achar;
                substitution[i]=achar;
                num++;
            }
        }
        return num;
}

//based on calibrated position-depended error rates
int SOLiDread::add_calib_error_1st(vector<short>&qual){
        int num=0;
        if(qual.size()<seq_read.size()) qual.resize(seq_read.size());
        for(size_t i=0; i<seq_read.size(); i++){
          int q=gsl_ran_poisson(gsl_p_1st[i],cal_qual_1st[i]);
          qual[i]=q;
          if(r_prob()<prob_err[q]){
                char achar=seq_read[i];
                while(seq_read[i]==achar){ achar=rand_base(); }
                seq_read[i]=achar;
                substitution[i]=achar;
                num++;
            }
        }
        return num;
}

int SOLiDread::add_calib_error_2nd(vector<short>&qual){
       	int num=0;
        if(qual.size()<seq_read.size()) qual.resize(seq_read.size());
        for(size_t i=0; i<seq_read.size(); i++){
          int q=gsl_ran_poisson(gsl_p_2nd[i],cal_qual_2nd[i]);
          qual[i]=q;
          if(r_prob()<prob_err[q]){
                char achar=seq_read[i];
                while(seq_read[i]==achar){ achar=rand_base(); }
                seq_read[i]=achar;
                substitution[i]=achar;
                num++;
            }
        }
        return num;
}


int SOLiDread::add_calib_error_1st_stat(vector<short>&err_pos){
        int num=0;
        for(size_t i=0; i<seq_read.size(); i++){
            if(r_prob()<cal_err_rate_1st[i]){
              err_pos[i]=1;
                char achar=seq_read[i];
                while(seq_read[i]==achar){ achar=rand_base(); }
                seq_read[i]=achar;
                substitution[i]=achar;
                num++;
            }
        }
        return num;
};

int SOLiDread::add_calib_error_2nd_stat(vector<short>&err_pos){
        int num=0;
        for(size_t i=0; i<seq_read.size(); i++){
            if(r_prob()<cal_err_rate_2nd[i]){
              err_pos[i]=1;
                char achar=seq_read[i];
                while(seq_read[i]==achar){ achar=rand_base(); }
                seq_read[i]=achar;
                substitution[i]=achar;
                num++;
            }
        }
        return num;
};


//based on empirical dist of quali scores
int SOLiDread::add_error(vector<short>&qual){
       	if(qual.size()!=seq_read.size()){
            cerr<<"error call in adding substitution\n";
            //cerr<<qual.size()<<"\t"<<seq_read.size()<<endl;
            return 0;
        }
        int num=0;
        for(size_t i=0; i<qual.size(); i++){
            if(r_prob()<prob_err[qual[i]]){
                char achar=seq_read[i];
                while(seq_read[i]==achar){ achar=rand_base(); }
                seq_read[i]=achar;
                substitution[i]=achar;
                num++;
            }
        }
        return num;
};

//based on empirical dist of quali scores
int SOLiDread::add_error_stat(vector<short>&qual, vector<long>&err_pos){
        if(qual.size()!=seq_read.size()){
            cerr<<"error call in adding substitution\n";
            //cerr<<qual.size()<<"\t"<<seq_read.size()<<endl;
            return 0;
        }
        int num=0;
        for(size_t i=0; i<qual.size(); i++){
            if(r_prob()<prob_err[qual[i]]){
              err_pos[i]+=1;
                char achar=seq_read[i];
                while(seq_read[i]==achar){ achar=rand_base(); }
                seq_read[i]=achar;
                substitution[i]=achar;
                num++;
            }
        }
        return num;
}

//base on bionomial error substitution rate
int SOLiDread::add_error(int read_len){
       	int sub_num=0;
        for(int i=(int)sub_rate.size()-1; i>=0; i--){
            if(sub_rate[i]>=r_prob()){
                sub_num=i+1;
                for(int j=i; j>=0;){
                    int pos=(int) floor(read_len*r_prob());
                    if(substitution.count(pos)==0){
                        char achar=seq_read[pos];
                        while(seq_read[pos]==achar){ achar=rand_base(); }
                        substitution[pos]=achar;
                        seq_read[pos]=achar;
                        j--;
                    }
                }
                break;
            }
        }
        return sub_num;
}


string SOLiDread::reverse_comp(){
       	string seq=seq_ref;
       	reverse(seq.begin(), seq.end());  
	for(int i=0; i<seq.length(); i++){
	       	if (seq[i] == 'a') seq[i] = 't';
	       	else if(seq[i] == 't') seq[i] = 'a';
	       	else if(seq[i] == 'c') seq[i] = 'g';
	       	else if(seq[i] == 'g') seq[i] = 'c';
	       	else seq[i] = 'n';
       	}
       	return seq;
}

