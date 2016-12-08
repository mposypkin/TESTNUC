all dep clean indent tests::
	cd testboxcon && $(MAKE) $@ && cd .. \\
	cd testdejong && $(MAKE) $@ && cd ..
	
doc: indent doxy

clean::
	rm -rf *~ PI* core bin/* obj/* tmp *.log
