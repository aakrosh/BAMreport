profile1dir=$(prefix)/Profiles/InDel_Profiles
profile1_DATA=Profiles/InDel_Profiles/*
profile2dir=$(prefix)/Profiles/Base-Calling_Profiles
profile2_DATA=Profiles/Base-Calling_Profiles/*
profile3dir=$(prefix)/Profiles/GC-depth_Profiles
profile3_DATA=Profiles/GC-depth_Profiles/*

AUTOMAKE_OPTIONS=gnu
bindir=$(prefix)/
bin_PROGRAMS=pirs
pirs_SOURCES=main.cpp simulate_Illumina_reads.cpp simulate_diploid_genome.cpp simulate.cpp load_file.cpp MaskQvalsByEamss.cpp gzstream.cpp
