
#include "Tools.h"

#include <fcntl.h>
#include <iomanip>
#include <unistd.h>
#include <cstring>
#include <stdio.h>

void Tools::printMemFootPrint(std::string tag) {

    int val[4];

    char filename[80];
    char sbuf[1024];
    char* S;
    int fd, num_read;
    long lmem;
    pid_t numpro;
    pid_t getpid(void);

    numpro = getpid();

    sprintf(filename, "/proc/%ld/status", (long)numpro);
    fd = open(filename, O_RDONLY, 0);
    num_read=read(fd,sbuf,(sizeof sbuf)-1);
    close(fd);

    // Peak resident set size
    S=strstr(sbuf,"VmRSS:")+8;
    val[1] = (int)atoi(S);
    // Peak virtual memory usage
    S=strstr(sbuf,"VmSize:")+8;
    val[2] = atoi(S);

    std::cout << "=== Mem usage === " <<  std::setw(20) << tag << "\t=== " << std::setw(6)
              << "\t VmRSS  << " <<  (int)((double)val[1]/1024.) << " Mo" << std::endl;

}

