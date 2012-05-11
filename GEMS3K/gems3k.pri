       HEADERS	 += $$GEMS3K_H/verror.h  \
                    $$GEMS3K_H/gdatastream.h  \
                    $$GEMS3K_H/v_user.h \
                    $$GEMS3K_H/tnt_i_refvec.h \
                    $$GEMS3K_H/tnt_array1d.h \
                    $$GEMS3K_H/tnt_array2d.h \
                    $$GEMS3K_H/tnt.h \
                    $$GEMS3K_H/jama_cholesky.h \
                    $$GEMS3K_H/jama_lu.h \
                    $$GEMS3K_H/num_methods.h \
                    $$GEMS3K_H/s_fgl.h \
#                    $$GEMS3K_H/s_sorption.h \
#                    $$GEMS3K_H/s_kinmet.h \
                    $$GEMS3K_H/m_param.h  \
                    $$GEMS3K_H/ms_multi.h \
                    $$GEMS3K_H/databr.h \
                    $$GEMS3K_H/datach.h \
                    $$GEMS3K_H/node.h \
                    $$GEMS3K_H/nodearray.h \
#                    $$GEMS3K_H/particlearray.h \
                    $$GEMS3K_H/io_arrays.h

        SOURCES	  +=  $$GEMS3K_CPP/gdatastream.cpp  \
                      $$GEMS3K_CPP/num_methods.cpp \
                      $$GEMS3K_CPP/s_fgl.cpp \
                      $$GEMS3K_CPP/s_fgl1.cpp \
                      $$GEMS3K_CPP/s_fgl2.cpp \
                      $$GEMS3K_CPP/s_fgl3.cpp \
                      $$GEMS3K_CPP/s_fgl4.cpp \
#                      $$GEMS3K_CPP/s_sorption.cpp \
#                      $$GEMS3K_CPP/s_kinmet.cpp \
                      $$GEMS3K_CPP/ipm_chemical.cpp \
                      $$GEMS3K_CPP/ipm_chemical2.cpp \
                      $$GEMS3K_CPP/ipm_chemical3.cpp \
                      $$GEMS3K_CPP/ipm_main.cpp \
                      $$GEMS3K_CPP/ipm_simplex.cpp \
                      $$GEMS3K_CPP/node.cpp \
                      $$GEMS3K_CPP/nodearray.cpp \
                      $$GEMS3K_CPP/node_format.cpp \
#                      $$GEMS3K_CPP/particlearray.cpp \
                      $$GEMS3K_CPP/ms_multi_file.cpp \
                      $$GEMS3K_CPP/ms_multi_format.cpp \
                      $$GEMS3K_CPP/io_arrays.cpp

