AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

AM_CPPFLAGS = -O3

bin_PROGRAMS = art_illumina art_454 art_SOLiD

art_illumina_SOURCES = art_illumina_src/art.h art_illumina_src/empdist.h art_illumina_src/readSeqFile.h art_illumina_src/seqRead.h art_illumina_src/samRead.h \
		       art_illumina_src/art_illumina.cpp art_illumina_src/art_qual_scale.cpp art_illumina_src/empdist.cpp art_illumina_src/readSeqFile.cpp art_illumina_src/seqRead.cpp art_illumina_src/samRead.cpp
art_454_SOURCES = art_454_src/art.h art_454_src/read_profile.h art_454_src/readSeqFile.h art_454_src/seqRead.h art_454_src/samRead.h \
		  art_454_src/art_454.cpp art_454_src/art.cpp art_454_src/read_profile.cpp art_454_src/readSeqFile.cpp art_454_src/seqRead.cpp  art_454_src/samRead.cpp

art_SOLiD_SOURCES= art_SOLiD_src/art.h art_SOLiD_src/readSeqFile.h art_SOLiD_src/read_SOLiD.h art_SOLiD_src/samRead.h art_SOLiD_src/utility.hpp \
	    art_SOLiD_src/art_SOLiD.cpp art_SOLiD_src/art.cpp art_SOLiD_src/readSeqFile.cpp art_SOLiD_src/read_SOLiD.cpp art_SOLiD_src/samRead.cpp

dist_bin_SCRIPTS = aln2bed.pl map2bed.pl ART_profiler_454/art_profiler_454 ART_profiler_illumina/combinedAvg.pl ART_profiler_illumina/fastqReadAvg.pl ART_profiler_illumina/art_profiler_illumina ART_profiler_illumina/empDist.pl ART_profiler_illumina/summation.pl
#dist_noinst_SCRIPTS = autogen.sh
