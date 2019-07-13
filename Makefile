# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.

include make.inc

# execute a target irrespective of the presence of a file or directory
# with the same name
.PHONY: ksum solver

default :
	@echo ''
	@echo 'To install WRAP, type at the shell prompt:'
	@echo 'make [-j] target'
	@echo ''
	@echo 'where target identifies one or multiple PACKAGES:'
	@echo '    ksum        DMFT outer loop package'
	@echo '  solver        CTQMC impurity solver'
	@echo '   clean        remove executables and objects'
	@echo '     all        compile all components'
	@echo ''

all  : ksum solver

ksum :
	@echo "Compiling ksum :"
	@if test -d ksum; then \
	(cd ksum ; $(MAKE) TLDEPS= ksum || exit 1); fi

solver :
	@echo "Compiling solver :"
	@if test -d solver; then \
	(cd solver ; $(MAKE) TLDEPS= all || exit 1); fi

# remove object files and executables
clean :
	@echo "Cleaning :"
	@touch make.inc
	@for dir in \
		ksum solver \
	; do \
	    if test -d $$dir ; then \
		( cd $$dir ; \
		$(MAKE) TLDEPS= clean ) \
	    fi \
	done