FC = g77

OBJ = daerfj.o \
    datrfj.o \
    dchqfj.o \
    dcpffj.o \
    dcprfj.o \
    dctsfj.o \
    dedffj.o \
    deptfg.o \
    depths.o \
    deptsp.o \
    dfdcfj.o \
    dfdcjs.o \
    dfdcsp.o \
    dficfj.o \
    dficjs.o \
    dficsp.o \
    dgdffj.o \
    dgl1fg.o \
    dgl1hs.o \
    dgl1sp.o \
    dgl2co.o \
    dgl2fg.o \
    dgl2hs.o \
    dgl2sp.o \
    dhhdfj.o \
    diacfj.o \
    diadfj.o \
    diaofj.o \
    diarfj.o \
    dierfj.o \
    dierjs.o \
    diersp.o \
    dljcfg.o \
    dminfg.o \
    dminhs.o \
    dminsp.o \
    dminxb.o \
    dmsabc.o \
    dmsafg.o \
    dmsahs.o \
    dmsasp.o \
    dodcfg.o \
    dodchs.o \
    dodcps.o \
    dodcsp.o \
    dpjbds.o \
    dpjbfg.o \
    dpjbhs.o \
    dpjbsp.o \
    dsfdfj.o \
    dsfdjs.o \
    dsfdsp.o \
    dsfifj.o \
    dsfijs.o \
    dsfisp.o \
    dsscfg.o \
    dsschs.o \
    dsscsp.o \
    dsysfv.o \
    dsysjs.o 

all: $(OBJ) lib

lib:
	gfortran -shared -o libminpack-2.so *.o
