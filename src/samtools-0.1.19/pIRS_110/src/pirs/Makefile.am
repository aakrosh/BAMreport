SUBDIRS = src
profile1dir=$(prefix)/Profiles/InDel_Profiles
profile1_DATA=$(SUBDIRS)/Profiles/InDel_Profiles/*
profile2dir=$(prefix)/Profiles/Base-Calling_Profiles
profile2_DATA=$(SUBDIRS)/Profiles/Base-Calling_Profiles/*
profile3dir=$(prefix)/Profiles/GC-depth_Profiles
profile3_DATA=$(SUBDIRS)/Profiles/GC-depth_Profiles/*

AUTOMAKE_OPTIONS=gnu
bindir=$(prefix)/
bin_PROGRAMS=pirs
pirs_SOURCES=$(SUBDIRS)/main.cpp $(SUBDIRS)/simulate_Illumina_reads.cpp $(SUBDIRS)/simulate_diploid_genome.cpp $(SUBDIRS)/simulate.cpp $(SUBDIRS)/load_file.cpp $(SUBDIRS)/MaskQvalsByEamss.cpp $(SUBDIRS)/gzstream.cpp
