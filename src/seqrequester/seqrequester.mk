TARGET   := seqrequester
SOURCES  := seqrequester.C \
            extract.C \
            generate.C \
            microsatellite.C \
            mutate.C \
            partition.C \
            sample.C \
            shiftregister.C \
            shiftregister-emit-fast.C \
            shiftregister-search-fast.C \
            shiftregister-search-slow.C \
            simulate.C \
            summarize.C

SRC_INCDIRS  := .

#  If we're part of Canu, build with canu support and use Canu's copy of
#  meryl-utility.  Otherwise, don't.  (meryl does this too)
ifneq ($(wildcard stores/sqStore.H), )
  SRC_CXXFLAGS := -DCANU
  SRC_INCDIRS  := ../../../utility/src/utility ../../../stores

#  If we're part of something else, include the something else's
#  utility directory.
else ifneq ($(wildcard seqrequester/src/seqrequester/seqrequester.C), )
  SRC_INCDIRS  := ../../../utility/src/utility

#  Otherwise, we're building directly in the seqrequester repo.
else
  SRC_INCDIRS  := ../utility/src/utility

endif

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -l${MODULE}
TGT_PREREQS := lib${MODULE}.a
