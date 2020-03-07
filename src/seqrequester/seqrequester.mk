
#  If 'make' isn't run from the root directory, we need to set these to
#  point to the upper level build directory.
ifeq "$(strip ${BUILD_DIR})" ""
  BUILD_DIR    := ../$(OSTYPE)-$(MACHINETYPE)/obj
endif
ifeq "$(strip ${TARGET_DIR})" ""
  TARGET_DIR   := ../$(OSTYPE)-$(MACHINETYPE)
endif

#            shiftregister-emit-fast.C \
#            shiftregister-emit-slow.C \

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
TGT_LDLIBS  := -lseqrequester
TGT_PREREQS := libseqrequester.a

SUBMAKEFILES :=
