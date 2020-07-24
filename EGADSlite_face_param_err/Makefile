LIBS = \
  -legadslite -L$(ESP_ROOT)/lib
  
IDIR = -I$(ESP_ROOT)/include

CC = gcc

all: egadslite_face_param_err

egadslite_face_param_err: egadslite_face_param_err.c
	$(CC) $? $(IDIR) $(LIBS) -o $@

clean:
	rm -f egadslite_face_param_err egadslite_face_param_err.o
