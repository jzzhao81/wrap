# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License. See the file `License' in the root directory
# of the present distribution.

# execute a target irrespective of the presence of a file or directory
# with the same name
.PHONY: base capi gardenia manjushaka

default: all

all: base capi gardenia manjushaka

base :
	@echo "Compiling solver base :"
	@if test -d base; then \
	(cd base ; $(MAKE) TLDEPS= all || exit 1); fi

capi :
	@echo "Compiling solver capi :"
	@if test -d capi; then \
	(cd capi ; $(MAKE) TLDEPS= all || exit 1); fi

gardenia :
	@echo "Compiling solver gardenia :"
	@if test -d gardenia; then \
	(cd gardenia ; $(MAKE) TLDEPS= all || exit 1); fi

manjushaka :
	@echo "Compiling solver manjushaka :"
	@if test -d manjushaka; then \
	(cd manjushaka ; $(MAKE) TLDEPS= all || exit 1); fi

# remove object files and executables
clean :
	@echo "Cleaning solver :"
	@for dir in \
		base capi gardenia manjushaka \
	; do \
	    if test -d $$dir ; then \
		( cd $$dir ; \
		$(MAKE) TLDEPS= clean ) \
	    fi \
	done