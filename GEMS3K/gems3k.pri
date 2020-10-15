SJSON_CPP = $$GEMS3K_H/../simdjson/singleheader
SJSON_H   = $$SJSON_CPP
DEPENDPATH += $$SJSON_H
INCLUDEPATH += $$SJSON_H

 HEADERS	 += $$GEMS3K_H/verror.h  \
                    $$GEMS3K_H/gdatastream.h  \
                    ##$$GEMS3K_H/v_user.h \
                    $$GEMS3K_H/tnt_i_refvec.h \
                    $$GEMS3K_H/tnt_array1d.h \
                    $$GEMS3K_H/tnt_array2d.h \
                    $$GEMS3K_H/tnt.h \
                    $$GEMS3K_H/jama_cholesky.h \
                    $$GEMS3K_H/jama_lu.h \
                    $$GEMS3K_H/num_methods.h \
                    $$GEMS3K_H/s_solmod.h \
                    $$GEMS3K_H/s_sorpmod.h \
                    $$GEMS3K_H/s_kinmet.h \
                    ##$$GEMS3K_H/m_param.h  \
                    $$GEMS3K_H/m_const_base.h  \
                    $$GEMS3K_H/ms_multi.h \
                    $$GEMS3K_H/databr.h \
                    $$GEMS3K_H/datach.h \
                    $$GEMS3K_H/node.h \
                    $$GEMS3K_H/nodearray.h \
#                    $$GEMS3K_H/particlearray.h \
                    #$$GEMS3K_H/io_arrays.h \
                    #$$GEMS3K_H/io_json.h \
                    $$GEMS3K_H/io_template.h \
                    $$GEMS3K_H/io_nlohmann.h \
                    $$GEMS3K_H/io_keyvalue.h \
                    $$GEMS3K_H/gems3k_impex.h \
                    $$GEMS3K_H/activities.h \
                    $$GEMS3K_H/kinetics.h \
                    $$GEMS3K_H/v_detail.h \
    $$PWD/io_simdjson.h \
                    $$SJSON_H/simdjson.h


        SOURCES	  +=  $$GEMS3K_CPP/gdatastream.cpp  \
                      $$GEMS3K_CPP/num_methods.cpp \
                      $$GEMS3K_CPP/s_solmod.cpp \
                      $$GEMS3K_CPP/s_solmod2.cpp \
                      $$GEMS3K_CPP/s_solmod3.cpp \
                      $$GEMS3K_CPP/s_solmod4.cpp \
                      $$GEMS3K_CPP/s_solmod5.cpp \
                      $$GEMS3K_CPP/s_solmod6.cpp \
                      $$GEMS3K_CPP/s_sorpmod.cpp \
                      $$GEMS3K_CPP/s_kinmet.cpp \
                      $$GEMS3K_CPP/ipm_chemical.cpp \
                      $$GEMS3K_CPP/ipm_chemical2.cpp \
                      $$GEMS3K_CPP/ipm_chemical3.cpp \
                      $$GEMS3K_CPP/ipm_chemical4.cpp \
                      $$GEMS3K_CPP/ipm_main.cpp \
                      $$GEMS3K_CPP/ipm_simplex.cpp \
                      $$GEMS3K_CPP/node.cpp \
                      $$GEMS3K_CPP/nodearray.cpp \
                      $$GEMS3K_CPP/nodearray_new.cpp \
                      $$GEMS3K_CPP/node_format.cpp \
#                      $$GEMS3K_CPP/particlearray.cpp \
                      $$GEMS3K_CPP/ms_multi_file.cpp \
                      $$GEMS3K_CPP/ms_multi_format.cpp \
                      #$$GEMS3K_CPP/io_arrays.cpp \
                      #$$GEMS3K_CPP/io_json.cpp \
                      $$GEMS3K_CPP/io_template.cpp \
                      $$GEMS3K_CPP/io_nlohmann.cpp \
                      $$GEMS3K_CPP/io_keyvalue.cpp \
                      $$GEMS3K_CPP/gems3k_impex.cpp \
                      $$GEMS3K_CPP/node_copy.cpp \
                      $$GEMS3K_CPP/ms_multi_copy.cpp \
                      $$GEMS3K_CPP/node_activities.cpp \
                      $$GEMS3K_CPP/node_kinetics.cpp \
                      $$GEMS3K_CPP/s_activity.cpp \
                      $$GEMS3K_CPP/s_activity2.cpp \
                      $$GEMS3K_CPP/s_activity3.cpp  \
                      #$$GEMS3K_CPP/ms_param.cpp \
                      $$GEMS3K_CPP/v_detail.cpp \
    $$PWD/io_simdjson.cpp \
                      $$SJSON_CPP/simdjson.cpp

