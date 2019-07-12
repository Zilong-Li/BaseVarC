bin_PROGRAMS = BaseVarC

BaseVarC_CPPFLAGS = \
		-I$(top_srcdir)/SeqLib \
		-I$(top_srcdir)/SeqLib/htslib \
		-I$(top_srcdir)/fmt/include

BaseVarC_LDADD = \
		$(top_builddir)/SeqLib/src/libseqlib.a \
		$(top_builddir)/SeqLib/htslib/libhts.a \
		$(LDFLAGS)

BaseVarC_SOURCES = BaseVarC.cpp BamProcess.cpp BaseType.cpp Algorithm.cpp
