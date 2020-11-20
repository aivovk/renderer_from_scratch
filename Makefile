CPP = g++ -Wall #-g -pg
CPPFLAGS     =
LIBS         = -lm -lSDL2

DESTDIR = ./
TARGET  = main

WORKDIR = ./build/

OBJECTS := $(patsubst %.cpp,$(WORKDIR)%.o,$(wildcard *.cpp))
DEPS := $(OBJECTS:%.o=%.d)

all: $(DESTDIR)$(TARGET)

$(DESTDIR)$(TARGET): $(OBJECTS)
	$(CPP) -o $(DESTDIR)$(TARGET) $(OBJECTS) $(LIBS)

$(WORKDIR)%.o: %.cpp
	$(CPP) $(CPPFLAGS) -MMD -MF $(WORKDIR)$*.d -c $(CFLAGS) $< -o $@

-include $(DEPS)

clean:
	-rm -f $(OBJECTS)
	-rm -f $(TARGET)
	-rm -f $(DEPS)
