/*! \mainpage SMILEI overview

 \section intro Introduction

 This is the introduction.

 You can add interesting cite commands like this \cite Grech2011 (you need to put the BIBTEX record in doc/smilei.bib)

 \section install_sec Installation

 \subsection Pre-requisites

 \subsection compiling

 This is (o should be) a massive parallel code which heavily uses MPI routines.
 The compilation works fine using <tt>g++</tt> and intel <tt>icpc</tt> compiler and openmpi.

 Two basic compilation mode: "release" and debug (accessible via the make commands <tt>make release</tt> and
 <tt>make debug</tt>).

 "debug" mode is without optimization and prints debug messages from <tt>DEBUG()</tt> macro (see \ref macros page).

 "release" mode is with <tt>-O3</tt> option and suppress the debug messages from <tt>DEBUG()</tt> macro.

 By default the code is compiled in "release" mode.

 \subsection Dependencies
 The code needs hdf5 libraries installed with parallel support (on mac with macport you need to install the package with:

 <tt>sudo port install hdf5-18 +openmpi</tt>

 */

/*! @file Pic.h

 @brief Pic.h

 @author tommaso vinci
 @date 2013-02-15
 */

#include <string>

#include "Tools.h"

//! verbosity level use verbose keyword in data parameter
#ifdef  __DEBUG
int debug_level = 10;
#else
int debug_level = 0;
#endif

//! main function
int main (int argc, char* argv[]);


void welcome() { // jp2a --width=80 -i smileiLogo.jpg | sed 's/^/MESSAGE(\"/'| sed 's/$/\");/' 
    MESSAGE("                                                                                ");
    MESSAGE("                        .','     .','.    .','.     ','.     .,,.     .,,.      ");
    MESSAGE("                      .dxxxxxc  oxxxxxl  lxxxxxo  cxxxxxd. ;xxxxxx. 'xxxxxx,    ");
    MESSAGE("                     .cxxxxxxx;;xxxxxxx:;xxxxxxx:,xxxxxxxc.xxxxxxxd.dxxxxxxx    ");
    MESSAGE("                     ,.dxxxxxc,.oxxxxxl'.lxxxxxo..cxxxxxd.';xxxxxx.,'xxxxxx,    ");
    MESSAGE("                 ....:..,cc:'.;..,:c:'.,..':c:,..,.':cc,..;..;cc;..;..;c;.      ");
    MESSAGE("              ,lddo; , 'ldxo:.' .cdxdc....cdxdc. ..:oxdl' ' ;oddl, ,            ");
    MESSAGE("             lxxxxxxd;;xxxxxxx;,xxxxxxx:,xxxxxxx;,xxxxxxx;,dxxxxxxl,            ");
    MESSAGE("             lxxxxxxd;:xxxxxxx:,xxxxxxx:,xxxxxxx;,xxxxxxx:,dxxxxxxl             ");
    MESSAGE("              ,oxxo: , ,lxxd:.' .ldxdc....cdxdl. ..:dxxl, ' :oxxo,              ");
    MESSAGE("               .,c;..:..,cc:..;........,'........,........;....                 ");
    MESSAGE("             .dxxxxx,,.dxxxxx:'        ..        .        '                     ");
    MESSAGE("             dxxxxxxx:cxxxxxxxc        ..        .        .                     ");
    MESSAGE("             'xxxxxx;,.dxxxxxc,        ..        .                              ");
    MESSAGE("               .;c:..:..;cc:'.;........,.........,....                          ");
    MESSAGE("              ,lddl; , 'cddo; ' .codo:....:odoc. . ;oddc'                       ");
    MESSAGE("             cxxxxxxo;;xxxxxxx;'xxxxxxx;,xxxxxxx;,xxxxxxx;                      ");
    MESSAGE("             lxxxxxxd;:xxxxxxx:,xxxxxxx:;xxxxxxx:,xxxxxxx:'                     ");
    MESSAGE("              ;oxxd:., ,oxxdc.' 'lxxxl....lxxxl' ..cdxxo, '                     ");
    MESSAGE("                 ....:..':c;..;..':c;..,...;c:'..,..;c:'..;..,:'.               ");
    MESSAGE("                     ,.o,...c:'.l;...:c'.c:...;l..:c...,o.',l...'o.             ");
    MESSAGE("                     .c,     xc::     ll:l     :c:x     ,c;x     .d             ");
    MESSAGE("                      .o.  .:c,.o'   ,o'.o,   'o..c:.  .o.':c.  .l'             ");
    MESSAGE("                        .,lc,.;..;llc,.,'.,cll;..,.,cll:..;.'cl;'               ");
    MESSAGE("                              '        ..        . ;:;::. ' ,:;::'              ");
    MESSAGE("                              '        ..        ,d.    c;,o'    :c             ");
    MESSAGE("                     ,        '        ..        ;x     ::,d.    ,l             ");
    MESSAGE("                     ,        '        ..        ..c:,;c; '.::;;::              ");
    MESSAGE("               .';,..:..';;,..;...;;;..,...;;;...,..,::'..;..,:'.               ");
    MESSAGE("             .l,..'l',.l;...c;' c:...c:..:c...:c..;c...;l.''l'..,l.             ");
    MESSAGE("            .o.     x;c,     dc::     lc:l     :c:d     ,c;x     .o             ");
    MESSAGE("            ',l.  .c:,.o.   ;l,.o'   ,o'.o,   'o.'l;   .o.':c.  .l,             ");
    MESSAGE("        ....;.'cloc'.:..:loc,.;..;lol;.,'.;lol;..,.,col:..;.'cl;'               ");
    MESSAGE("     ,::::' ' '::::, , .::::; ' .::::;....;::::. . ;::::. '                     ");
    MESSAGE("    o,    :c'c:    ,o;,l    .d;'o    .d;'d.    o,,d.    l,'                     ");
    MESSAGE("    d.    'o o'    .d.::     x.,l     o''o     l,.x     ::                      ");
    MESSAGE("    .c:,,::   ::,,:c.  ;c;,:c.  ,c;,;c'  'c;,;c,  .c:,;c;                       ");
    MESSAGE("                                                                                ");
    MESSAGE("         .;ccc;.                          ;lc.                      ..          ");
    MESSAGE("        ldddddddo                         ddd.                     dxxc         ");
    MESSAGE("       ;ddd.  ',.  .'. ..''   .''.   lxx, ddd.    .,;,'.   'xxo     lxxl        ");
    MESSAGE("       .oddd:'.    dddodddddldddddc  ,lc. ddd.  ;dddodddl. .cl;      dxx'       ");
    MESSAGE("         'cddddo;  ddd,  cddo  .ddd. .;;  ddd. ;ddo. .cddl  ;:'      ;xxc       ");
    MESSAGE("            .,oddl ddd.  ,ddl   ddd. lxx; ddd. odddoooooo: .xxx      ;xxc       ");
    MESSAGE("        ''.  .lddo ddd.  ,ddl   ddd. lxx; ddd. :ddc    .'. .xxx      dxx'       ");
    MESSAGE("       ;ddddddddo. ddd.  ,ddl   ddd. 'lc. ddd.  cdddolddd:  cl;     lxxl        ");
    MESSAGE("        .,:cc:;.   ';,    ;;.   .;,       .;,    .,:cc:,.          dxxc         ");
    MESSAGE("                                                                    ..          ");
    MESSAGE("                                                                                ");
}