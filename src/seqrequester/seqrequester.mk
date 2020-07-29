TARGET   := seqrequester
SOURCES  := seqrequester.C \
            extract.C \
            generate.C \
            mutate.C \
            sample.C \
            shiftregister.C \
            shiftregister-emit-fast.C \
            shiftregister-search-fast.C \
            shiftregister-search-slow.C \
            simulate.C \
            summarize.C

SRC_INCDIRS  := . .. ../utility/src/utility

TGT_LDFLAGS := -L${TARGET_DIR}/lib
TGT_LDLIBS  := -l${MODULE}
TGT_PREREQS := lib${MODULE}.a
