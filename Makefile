AR       = ar
ARFLAGS  = rus

CXX      = g++

CXXFLAGS = -g -MMD -MP -Wall -Wextra -Winit-self -Wno-unused-parameter -std=c++11 -O3

RM       = rm -f
LDFLAGS  = 
LIBS     =
INCLUDE  = -I../include/boost_1_69_0/
TARGET   = ./train
OBJDIR   = ./obj
SOURCES  = $(wildcard *.cc)
OBJECTS  = $(addprefix $(OBJDIR)/, $(SOURCES:.cc=.o))
DEPENDS  = $(OBJECTS:.o=.d)
$(TARGET): $(OBJECTS) $(LIBS)
	$(CXX) -o $@ $^ $(LDFLAGS)

#compile
$(OBJDIR)/%.o: %.cc
	@if [ ! -d $(OBJDIR) ];\
	then echo "mkdir -p $(OBJDIR)";mkdir -p $(OBJDIR);\
	fi
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

#clean and build
.PHONY: all
all: clean $(TARGET)

#clean
.PHONY:clean
clean:
	$(RM) $(OBJECTS) $(DEPENDS) $(TARGET)

-include $(DEPENDS)
