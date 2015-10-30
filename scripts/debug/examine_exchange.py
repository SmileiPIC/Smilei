import scipy

def treatline(b, patchcountarray):
    b = b.split()
    if(len(b) > 0):
        if (b[0][:11] == "patch_count") and (b[3] == "="):
            patchcountarray.append(float(b[4]))

def treatline2(b, send_mat, recv_mat, create_mat):
    b = b.split()
    if(len(b) > 1):
        if (b[1] == "creating") and (b[6] == "="):
            create_mat[int(b[0])] += 1
    if(len(b) > 2):
        if (b[2] == "sending"):
            send_mat[int(b[1]),int(b[-1])] += 1
    if(len(b) > 2):
        if (b[2] == "receiving"):
            recv_mat[int(b[1]),int(b[-1])] += 1
    
directories = [""]
logfiles = ["exchange_data.txt"]

for ifile in range(len(directories)):

    patchcountarray = []
    datafile = open(directories[ifile]+logfiles[ifile])
    
    location = 0
    location_old = -1
    nline = 0
    
    while (location != location_old):
        location_old = location
        buffer = datafile.readline()
        nline = nline + 1
        location = datafile.tell()
        treatline(buffer, patchcountarray)
    
    datafile.close()
    
    print nline, "lines were read, ", len(patchcountarray), " times were collected."
    
    patchcountarray = scipy.array(patchcountarray)


print patchcountarray 

nmpi = len(patchcountarray)

send_mat=scipy.zeros((nmpi,nmpi))
recv_mat=scipy.zeros((nmpi,nmpi))
create_mat=scipy.zeros((nmpi))

for ifile in range(len(directories)):

    datafile = open(directories[ifile]+logfiles[ifile])
    
    location = 0
    location_old = -1
    nline = 0
    
    while (location != location_old):
        location_old = location
        buffer = datafile.readline()
        nline = nline + 1
        location = datafile.tell()
        treatline2(buffer, send_mat, recv_mat, create_mat)
    
    datafile.close()

print create_mat

for i in range(nmpi):
    senttoi = send_mat[:,i].sum()
    if senttoi != create_mat[i]:
        print "error, ", sentoi, " patches are sent to ",i, "only ", create_mat[i]," patches are created"
    for j in range(nmpi):
        if send_mat[i,j] != recv_mat[j,i]:
            print "error, ",i, " sends ",send_mat[i,j], "to ",j
            print j, "receives ", recv_mat[j,i], "from ",i 

