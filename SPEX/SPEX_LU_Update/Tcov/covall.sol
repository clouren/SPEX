#!/bin/csh
        tcov -x cm.profile SPEX_*.c spex_*.c >& /dev/null
        echo -n "statments not yet tested: "
        ./covs > covs.out
        grep "#####" *tcov | wc -l
