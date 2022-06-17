CPPFLAGS=`root-config --cflags` 
LDFLAGS=-g -L$(ROOTSYS)/lib
LDLIBS=-lRooFit -lHtml -lMinuit -lRooFitCore `root-config --glibs` 
CC=g++

TARGET= app2

#all: $(TARGET)
all: $(TARGET) $(TARGET2)

$(TARGET): $(TARGET).cc
	$(CC) $(LDFLAGS) -o $(TARGET) $(TARGET).cc dirFit.cc wbEvent.cc wbPDF.cc $(CPPFLAGS) $(LDLIBS)
clean:
	$(RM) $(TARGET) $(TARGET2)
