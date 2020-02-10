LIBS = \
  -legads -L$(ESP_ROOT)/lib

CC = gcc

all: convert

convert: egadsExport.c
	$(CC) $? $(LIBS) -o $@

clean:
	rm -f convert egadsExport.o
