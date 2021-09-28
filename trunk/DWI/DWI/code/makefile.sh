FFLAGS = -O1


tracemodel:	 tracemain_table.o trace_time.o xkintgrl.o source.o fourt.o bessl_table3.o tracegreen_table.o dmatrx.o gmatrx.o ematrx.o xeqy.o cmatmult.o analinv.o 
	f77 $(FFLAGS) -o tracemodel tracemain_table.o trace_time.o xkintgrl.o source.o fourt.o bessl_table3.o tracegreen_table.o dmatrx.o gmatrx.o ematrx.o xeqy.o cmatmult.o analinv.o 

bholemod:	 tracemain_table.class.o trace_time.o xkintgrl.o source.o fourt.o bessl_table3.o tracegreen_table.o dmatrx.o gmatrx.o ematrx.o xeqy.o cmatmult.o analinv.o 
	f77 $(FFLAGS) -o bholemod tracemain_table.class.o trace_time.o xkintgrl.o source.o fourt.o bessl_table3.o tracegreen_table.o dmatrx.o gmatrx.o ematrx.o xeqy.o cmatmult.o analinv.o 

bessrun:	bessrun.o bfcns.o
	f77 $(FFLAGS) -o bessrun bessrun.o bfcns.o
clean:
	rm *.o

