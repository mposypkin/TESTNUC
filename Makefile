all dep clean indent tests::
	cd testboxcon && $(MAKE) $@ && cd .. \\
	cd testdejong && $(MAKE) $@ && cd .. \\
	cd testackley1 && $(MAKE) $@ && cd ..
	
doc: indent doxy

clean::
	rm -rf *~ PI* core bin/* obj/* tmp *.log
