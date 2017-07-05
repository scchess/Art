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

#pragma once

#include <cmath>
#include <ctime>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iterator>
#include <map>
#include <cstdio>
#include <set>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
//#include <gsl/gsl_sys.h>
//#include <gsl/gsl_machine.h>
//#include <gsl/gsl_precision.h>
//#include <gsl/gsl_nan.h>
//#include <gsl/gsl_pow_int.h>
#include <gsl/gsl_cdf.h>

//#include "empdist.h"
#include "utility.hpp"

using namespace std;
#define max_qual_value 80

class SOLiDread{
public:
    static char rand_base();
    static double prob_err[max_qual_value];
    static void set_err_prob(){
        for(int i=0; i<max_qual_value; i++){
            prob_err[i]=pow(10,-i/(double)10);
        }
    };
 
    vector <double> cal_qual_1st;
    vector <double> cal_qual_2nd;
    vector<gsl_rng*> gsl_p_1st;
    vector<gsl_rng*> gsl_p_2nd;

    vector <double> cal_err_rate_1st;
    vector <double> cal_err_rate_2nd;
    vector <double> ins_rate; //Bionomial cum_prob gsl_cdf_binomial_Q(unsigned int k, double p, unsigned int n) 
    vector <double> del_rate; //Binomial
    vector<double> sub_rate; //Binomial

    //static bool with_indel;
    int get_indel(int read_len);
    map<int,char,less<int> > indel;
    map<int,char> substitution;
    bool is_plus_strand;
    unsigned long bpos; //parent
    string seq_read;
    string seq_ref;

    map<char, map<char, char> > cs2seq;
    map<char, map<char, char> > seq2cs;
    vector< map <char, base_quality> > error_profile;


//    SOLiDread();
    SOLiDread(string profile);

    bool read_error_profile(string infile_profile);

    //read position and color dependent error profile for paired-end reads
    bool read_error_profile(istream& input);
    void ini_ran_qual();
    void ini_ran_qual(unsigned int r_seed);
    void set_rate(int read_len, double p, int max_num, vector <double>& rate);

    //convert base-space to color space, and incorporate sequencing errors
    void convert_seq2cs(string& c_space, vector<short>& qual, map<int,char>& error_pos, bool is_F3=true);

    //convert base-space to color space
    void convert_seq2cs(string& b_space, string& c_space);

    //convert color space to base space
//    void convert_cs2seq(string& c_space, string& b_space);
    void convert_cs2seq(string& c_space,  string& b_space, char base1st);

    //string aln_read;
    //string aln_ref;
    bool get_aln(string& aln_read, string& aln_ref);

    void ref2read();
    
    //based on based on calibrated position-depended error rates
    int add_calib_error_1st();

    int add_calib_error_2nd();

    //based on calibrated position-depended error rates
    int add_calib_error_1st(vector<short>&qual);
    int add_calib_error_2nd(vector<short>&qual);

    int add_calib_error_1st_stat(vector<short>&err_pos);
    int add_calib_error_2nd_stat(vector<short>&err_pos);

    //based on empirical dist of quali scores
    int add_error(vector<short>&qual);

    //based on empirical dist of quali scores
    int add_error_stat(vector<short>&qual, vector<long>&err_pos);

    //base on bionomial error substitution rate
    int add_error(int read_len);
    void clear(){
        indel.clear();
        substitution.clear();
        seq_read.clear();
        seq_ref.clear();
    };

    string reverse_comp();

};
