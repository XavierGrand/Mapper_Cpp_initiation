CXX=g++
CXXFLAGS=-Wall -ansi -pedantic
DEPEND=.deps
PROG=mapper
SOURCES= Common.cpp fastxfile.cpp fastxseq.cpp EncodedSequence.cpp SA_ES.cpp Main.cpp
HEADERS= Common.h fastxfile.h fastxseq.h EncodedSequence.h SA_ES.h
OBJECTS=$(SOURCES:.cpp=.o)
.PHONY: all clean
.SUFFIXES: .cpp .o .h

ifeq ((wildcard $(DEPEND)), )
	all: $(DEPEND)
		make
else
	all: $(PROG)
	include $(DEPEND)
endif

$(DEPEND): $(SOURCES)
	$(CXX) -MM $^ > $@

$(PROG): $(OBJECTS)
	$(CXX) $^ -o $@
.cpp.o:
	$(CXX) $(CXXFLAGS) $< -c

clean: 
	rm -f $(PROG)
	rm -f $(OBJECTS)
	rm -f $(DEPEND)

archive: $(HEADERS) $(SOURCES) Makefile
	tar -czf $(PROG).tgz $^
