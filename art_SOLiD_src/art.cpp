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

#include "art.h"
using namespace std;
bool art::next_read_indel(SOLiDread& a_read){
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
//amplicon reads
bool art::amp_read(SOLiDread& a_read){
    long pos=(long)0; // floor(r_prob()*valid_region); //pos in [0 ..len-1]           
    //string seq_ref;
    a_read.is_plus_strand=true;
//    if(r_prob()>0.5){
//        a_read.is_plus_strand=false;
//    }
//    if(a_read.is_plus_strand){
        a_read.seq_ref=ref_seq.substr(pos, read_len);
//    }
//    else{
//        a_read.seq_ref=ref_seq_cmp.substr(pos, read_len);
//    }
    a_read.bpos=pos;
    a_read.seq_read=a_read.seq_ref;
    return true;
}

//amplicon matepair-end read: the 1st and 2nd reads are the same strand 
bool art::amp_mate_read(SOLiDread& read_1, SOLiDread& read_2){
    long pos_1=(long) 0; 
    long pos_2=ref_seq.length()-read_len;
    bool is_plus_strand=true;
    read_1.seq_ref=ref_seq.substr(pos_1, read_len);
    read_2.seq_ref=ref_seq.substr(pos_2, read_len);
    read_1.is_plus_strand=is_plus_strand;
    read_1.bpos=pos_1;
    read_2.is_plus_strand=is_plus_strand;
    read_2.bpos=pos_2;
    return true;
}

//amplicon pair-end read: the second read is the reverse complemenaty strand 
bool art::amp_PE_read(SOLiDread& read_1, SOLiDread& read_2){
    long pos_1=(long)0;
    long pos_2=0;
    bool is_plus_strand=true;
    read_1.is_plus_strand=false;
    read_1.seq_ref=ref_seq_cmp.substr(pos_1, read_len_F5); 
    read_2.is_plus_strand=true;
    read_2.seq_ref=ref_seq.substr(pos_2, read_len);
    read_1.bpos=pos_1;
    read_2.bpos=pos_2;
    return true;
}

bool art::next_read(SOLiDread& a_read){
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

//matepair-end read: the 1st and 2nd reads are the same strand 
bool art::next_pair_read(SOLiDread& read_1, SOLiDread& read_2){
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

//pair-end read: the second read is the reverse complemenaty strand 
bool art::next_PE_read(SOLiDread& read_1, SOLiDread& read_2){
    int fragment_len=gaussain_mean+ (int)floor(gsl_ran_gaussian(gsl_R, gaussain_sigma));
    while (fragment_len<read_len || fragment_len>ref_seq.length()){
        fragment_len=gaussain_mean+ (int)floor(gsl_ran_gaussian(gsl_R, gaussain_sigma));
    }
    long pos_1=(long) floor((ref_seq.length()-fragment_len)*r_prob())+fragment_len-read_len_F5;
    long pos_2=ref_seq.length()-(pos_1+read_len+read_len_F5-fragment_len);
    bool is_plus_strand=true;
    if(r_prob()>0.5){
        is_plus_strand=false;
    }
    if(is_plus_strand){
        read_1.is_plus_strand=false;
        read_1.seq_ref=ref_seq_cmp.substr(pos_1, read_len_F5); 
        read_2.is_plus_strand=true;
        read_2.seq_ref=ref_seq.substr(pos_2, read_len);
    }
    else{
        read_1.is_plus_strand=true;
        read_1.seq_ref=ref_seq.substr(pos_1, read_len_F5);
        read_2.is_plus_strand=false;
        read_2.seq_ref=ref_seq_cmp.substr(pos_2, read_len);
    }
    read_1.bpos=pos_1;
    read_2.bpos=pos_2;
    return true;
}

bool art::next_pair_read_indel(SOLiDread& read_1, SOLiDread& read_2){
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
bool art::next_pair_read_indel_cmp(SOLiDread& read_1, SOLiDread& read_2){
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
