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
#include "read_profile.h"

gsl_rng* art::gsl_R;
int art::gaussain_mean;
double art::gaussain_sigma;
double read_profile::prob_err[HIGHEST_QUAL];

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
    if(read_len>ref_seq.size()){
	    read_len=ref_seq.size()-10;
	    if(read_len<0) read_len=ref_seq.size();
	    valid_region=ref_seq.size()-read_len;
    }
   
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
    //a_read.seq_read=a_read.seq_ref;
    return true;
}

//454 paired-end reads
//bool art::next_pair_read(seqRead& read_1, seqRead& read_2, int len_of_read){
/*
bool art::next_pair_read_old(seqRead& read_1, seqRead& read_2){
    if(read_len <2*min_read_len) read_len=2*min_read_len; 

    int fragment_len=gaussain_mean+ (int)floor(gsl_ran_gaussian(gsl_R, gaussain_sigma));
    if (fragment_len<read_len || fragment_len>ref_seq.length()) return false;
    int len_1st=(int) floor((read_len-2*min_read_len)*r_prob())+min_read_len;
    int len_2nd=read_len-len_1st;
    long pos_1=(long) floor((ref_seq.length()-fragment_len)*r_prob());
    long pos_2=pos_1+fragment_len-len_2nd;
    bool is_plus_strand=true;
    if(r_prob()>0.5){
        is_plus_strand=false;
    }
    if(is_plus_strand){
        read_1.seq_ref=ref_seq.substr(pos_1, len_1st);
        read_2.seq_ref=ref_seq.substr(pos_2, len_2nd);
    }
    else{
        read_1.seq_ref=ref_seq_cmp.substr(pos_1, len_1st);
        read_2.seq_ref=ref_seq_cmp.substr(pos_2, len_2nd);
    }
    read_1.is_plus_strand=is_plus_strand;
    read_1.bpos=pos_1;
    read_2.is_plus_strand=is_plus_strand;
    read_2.bpos=pos_2;
    //cout<<pos_1<<" a: "<<read_1.seq_read<<endl<<pos_2<<" b: "<<read_2.seq_read<<endl;
    return true;
}
*/

//454 paired-end reads
bool art::next_pair_read(seqRead& read_1, seqRead& read_2){
//  read_len=len_of_read;
    //make sure each read has the minimum length 
    if(read_len <2*min_read_len) read_len=2*min_read_len; 

    int fragment_len=gaussain_mean+ (int)floor(gsl_ran_gaussian(gsl_R, gaussain_sigma));
    if (fragment_len<read_len || fragment_len>ref_seq.length()) return false;
/*    while (fragment_len<read_len || fragment_len>ref_seq.length()){
        fragment_len=gaussain_mean+ (int)floor(gsl_ran_gaussian(gsl_R, gaussain_sigma));
    }
 */
    int len_1st=(int) floor((read_len-2*min_read_len)*r_prob())+min_read_len;
    int len_2nd=read_len-len_1st;
    long pos_1=(long) floor((ref_seq.length()-fragment_len)*r_prob());
    long pos_2=ref_seq.length()-pos_1-fragment_len;
    bool is_plus_strand=true;
    if(r_prob()>0.5){
        is_plus_strand=false;
    }
    if(is_plus_strand){
        read_1.seq_ref=ref_seq.substr(pos_1, len_1st);
       	read_1.is_plus_strand=true;
        read_2.seq_ref=ref_seq_cmp.substr(pos_2, len_2nd);
       	read_2.is_plus_strand=false;
    }
    else{
        read_1.seq_ref=ref_seq_cmp.substr(pos_1, len_1st);
        read_1.is_plus_strand=false;
        read_2.seq_ref=ref_seq.substr(pos_2, len_2nd);
        read_2.is_plus_strand=true;
    }
    read_1.bpos=pos_1;
    read_2.bpos=pos_2;
    //cout<<pos_1<<" a: "<<read_1.seq_read<<endl<<pos_2<<" b: "<<read_2.seq_read<<endl;
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

//amplicon reads
bool art::amp_read(seqRead& a_read){
    if(read_len>ref_seq.size()){
	    read_len=ref_seq.size()-10;
	    if(read_len<0) read_len=ref_seq.size();
	    valid_region=ref_seq.size()-read_len;
    }
   
    long pos=(long) 0; 
    //string seq_ref;
    a_read.is_plus_strand=true;
        a_read.seq_ref=ref_seq.substr(pos, read_len);
    a_read.bpos=pos;
    //a_read.seq_read=a_read.seq_ref;
    return true;
}

//454 amplicon paired-end reads
//bool art::next_pair_read(seqRead& read_1, seqRead& read_2, int len_of_read){
bool art::amp_pair_read(seqRead& read_1, seqRead& read_2){
//  read_len=len_of_read;
    //make sure each read has the minimum length 
    if(read_len <2*min_read_len) read_len=2*min_read_len; 
    if (read_len>ref_seq.length()){
	    read_len=ref_seq.length();
	    if(read_len <2*min_read_len){
		    return false;
	    }
    }

    int len_1st=(int) floor((read_len-2*min_read_len)*r_prob())+min_read_len;
    int len_2nd=read_len-len_1st;
    long pos_1=(long) 0;
    long pos_2=(long) 0; 
    bool is_plus_strand=true;
    if(is_plus_strand){
        read_1.seq_ref=ref_seq.substr(pos_1, len_1st);
       	read_1.is_plus_strand=true;
        read_2.seq_ref=ref_seq_cmp.substr(pos_2, len_2nd);
       	read_2.is_plus_strand=false;
    }
    read_1.bpos=pos_1;
    read_2.bpos=pos_2;
    //cout<<pos_1<<" a: "<<read_1.seq_read<<endl<<pos_2<<" b: "<<read_2.seq_read<<endl;
    return true;
}
