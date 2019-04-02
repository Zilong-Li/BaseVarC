bin_PROGRAMS = BaseVarC

test_CPPFLAGS = \
		-I$(top_srcdir)/SeqLib \
		-I$(top_srcdir)/SeqLib/htslib 

test_LDADD = \
		$(top_builddir)/SeqLib/src/libseqlib.a \
		$(top_builddir)/SeqLib/htslib/libhts.a \
		$(LDFLAGS)

test_SOURCES = BaseVarC.cpp BamProcess.cpp
