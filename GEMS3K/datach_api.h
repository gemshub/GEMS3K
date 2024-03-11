//--------------------------------------------------------------------
// $Id$
//
/// \file datach_api.h
/// Interface for writing/reading DBR and DCH I/O files of GEMS3K
/// Functions that maintain DATACH and DATABR memory allocation
//
// Copyright (c) 2023 S.Dmytriyeva, D.Kulik
// <GEMS Development Team, mailto:gems2.support@psi.ch>
//
// This file is part of the GEMS3K code for thermodynamic modelling
// by Gibbs energy minimization <http://gems.web.psi.ch/GEMS3K/>
//
// GEMS3K is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation, either version 3 of
// the License, or (at your option) any later version.

// GEMS3K is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with GEMS3K code. If not, see <http://www.gnu.org/licenses/>.
//-------------------------------------------------------------------


#ifndef DATACH_API_H
#define DATACH_API_H

#include "datach.h"
#include "databr.h"
#include "gems3k_impex.h"

class GemDataStream;

namespace  dbr_dch_api {

/// Writes CSD (DATACH structure) to a text DCH file
/// \param brief_mode - Do not write data items that contain only default values
/// \param with_comments - Write files with comments for all data entries or as "pretty JSON"
template<typename TIO>
void datach_to_text_file(const DATACH* pCSD, TIO& out_format, bool use_thermofun, bool with_comments = true, bool brief_mode = false );

/// Reads CSD (DATACH structure) from a text DCH file
template<typename TIO>
void datach_from_text_file(DATACH* pCSD, TIO& in_format, bool use_thermofun);

/// Writes work node (DATABR structure) to a text DBR file
/// \param pCSD - Main chemical system data structure CSD
/// \param brief_mode - Do not write data items that contain only default values
/// \param with_comments - Write files with comments for all data entries or as "pretty JSON"
template<typename TIO>
void databr_to_text_file(const DATACH* pCSD, const DATABR* pCNode, TIO& out_format, bool with_comments = true, bool brief_mode = false);

/// Reads work node (DATABR structure) from a text DBR file
/// \param pCSD - Main chemical system data structure CSD
template<typename TIO>
void databr_from_text_file(const DATACH* pCSD, DATABR* pCNode, TIO& in_format);


/// Writes CSD (DATACH structure) to a json/key-value string
/// \param set_name - The same GEMS3K input set (currently GEM-Selektor asks for the "set" name as .lst file name)
/// \param brief_mode - Do not write data items that contain only default values
/// \param with_comments - Write files with comments for all data entries or as "pretty JSON"
std::string datach_to_string(const std::string& set_name, const DATACH* pCSD, bool with_comments, bool brief_mode);

/// Reads CSD (DATACH structure) from a json/key-value string
/// \param set_name - The same GEMS3K input set (currently GEM-Selektor asks for the "set" name as .lst file name)
bool datach_from_string(const std::string& set_name, DATACH* pCSD, const std::string data);

/// Writes work node (DATABR structure) to a json/key-value string
/// \param set_name - The same GEMS3K input set
/// \param pCSD - Main chemical system data structure CSD
/// \param brief_mode - Do not write data items that contain only default values
/// \param with_comments - Write files with comments for all data entries or as "pretty JSON"
std::string databr_to_string(const std::string& set_name, const DATACH* pCSD, const DATABR* pCNode, bool with_comments, bool brief_mode);

/// Reads work node (DATABR structure) from a json/key-value string
/// \param set_name - The same GEMS3K input set
/// \param pCSD - Main chemical system data structure CSD
bool databr_from_string(const std::string& set_name, const DATACH* pCSD, DATABR* pCNode, const std::string data);


/// Writes the contents of the work instance of the DATACH structure into a stream.
///   \param set_name - The same GEMS3K input set
///   \param stream    string or file stream.
///   \param type_f    defines if the file is in binary format (1), in text format (0) or in json format (2).
///   \param with_comments (text format only): defines the mode of output of comments written before each data tag and  content
///                 in the DBR file. If set to true (1), the comments will be written for all data entries (default).
///                 If   false (0), comments will not be written;
///                         (json format): interpret the flag with_comments=on as "pretty JSON" and
///                                   with_comments=off as "condensed JSON"
///  \param brief_mode     if true, tells that do not write data items,  that contain only default values in text format
void write_dch_format_stream(const std::string& set_name, const DATACH* pCSD, std::iostream& stream,
                             GEMS3KGenerator::IOModes type_f, bool with_comments, bool brief_mode);

///  Reads the contents of the work instance of the DATACH structure from a stream.
///   \param set_name - The same GEMS3K input set
///   \param stream    string or file stream.
///   \param type_f    defines if the file is in binary format (1), in text format (0) or in json format (2).
void read_dch_format_stream(const std::string& set_name, DATACH* pCSD, std::iostream& stream, GEMS3KGenerator::IOModes  type_f);

/// Writes the contents of the work instance of the DATABR structure into a stream.
///   \param set_name - The same GEMS3K input set
///   \param pCSD - Main chemical system data structure CSD
///   \param stream    string or file stream.
///   \param type_f    defines if the file is in binary format (1), in text format (0) or in json format (2).
///   \param with_comments (text format only): defines the mode of output of comments written before each data tag and  content
///                 in the DBR file. If set to true (1), the comments will be written for all data entries (default).
///                 If   false (0), comments will not be written;
///                         (json format): interpret the flag with_comments=on as "pretty JSON" and
///                                   with_comments=off as "condensed JSON"
///  \param brief_mode     if true, tells that do not write data items,  that contain only default values in text format
void write_dbr_format_stream(const std::string& set_name, const DATACH* pCSD, const DATABR* pCNode, std::iostream& stream,
                              GEMS3KGenerator::IOModes type_f,  bool with_comments, bool brief_mode);

///  Reads the contents of the work instance of the DATABR structure from a stream.
///   \param set_name - The same GEMS3K input set
///   \param pCSD - Main chemical system data structure CSD
///   \param stream    string or file stream.
///   \param type_f    defines if the file is in binary format (1), in text format (0) or in json format (2).
void read_dbr_format_stream(const std::string& set_name, const DATACH* pCSD, DATABR* pCNode, std::iostream& stream, GEMS3KGenerator::IOModes type_f);


/// Writes DATACH structure to a binary DCH file.
void datach_to_file(const DATACH* CSD, GemDataStream& ff);
/// Reads DATACH structure from a binary DCH file.
void datach_from_file(DATACH* CSD, GemDataStream& ff);
/// Writes work DATABR structure to a binary DBR file.
void databr_to_file(const DATACH* CSD, const DATABR* CNode, GemDataStream& ff);
/// Reads work DATABR structure from a binary DBR file.
void databr_from_file(const DATACH* CSD, DATABR* CNode, GemDataStream& ff);


/// Returns number of temperature and pressure grid points for one dependent component
long int gridTP(const DATACH* pCSD);
/// Checks if given temperature TK and pressure P fit within the interpolation
/// intervals of the DATACH lookup arrays (returns empty message) or not (returns error message)
std::string check_TP(const DATACH* CSD, double TK, double P);
/// Tests TK as a grid point for the interpolation of thermodynamic data.
/// \return index in the lookup grid array or -1  if it is not a grid point
long int check_grid_T(const DATACH* CSD, double TK);
/// Tests P as a grid point for the interpolation of thermodynamic data.
/// \return index in the lookup grid array or -1 if it is not a grid point
long int check_grid_P(const DATACH* CSD, double P);
/// Tests TK (K) and P (Pa) as a grid point for the interpolation of thermodynamic data using DATACH
/// lookup arrays. \return -1L if interpolation is needed, or 1D index of the lookup array element
/// if TK and P fit within the respective tolerances.
long int check_grid_TP(const DATACH* CSD, double TK, double P) ;

// Functions that maintain DATACH and DATABR memory allocation

/// Set default values(zeros) for DATACH structure
void datach_reset(DATACH* CSD);
/// Allocating memory for DATACH structure
void datach_realloc(DATACH* pCSD);
/// Freeing memory for DATACH structure
void datach_free(DATACH* pCSD);

/// Set default values(zeros) for DATABR structure
void databr_reset(DATABR *pCNode, long int level);
/// Allocating memory for DATABR structure
///   \param pCSD - dimensions nICb, nDCb, nPHb, nPSb see in the DATACH structure
void databr_realloc(const DATACH* pCSD, DATABR* pCNode);
/// Freeing memory for DATABR structure
void databr_free_internal(DATABR *pCNode);

}  // dbr_dch_api

#endif // DATACH_API_H
