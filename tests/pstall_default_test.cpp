#include "datach_api.h"
#include "ms_multi.h"

#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>

namespace
{
std::string read_file(const std::string& path)
{
    std::ifstream input(path);
    if (!input)
        throw std::runtime_error("Cannot open test input: " + path);

    return {std::istreambuf_iterator<char>(input), std::istreambuf_iterator<char>()};
}
}

int main()
{
    DATACH dch{};
    dbr_dch_api::datach_reset(&dch);

    try
    {
        const std::string data_dir = GEMS3K_TEST_DATA_DIR;
        std::fstream dch_input(data_dir + "/CalcColumn-dch.dat", std::ios::in);
        dbr_dch_api::read_dch_format_stream(
            "pstall_default_test", &dch, dch_input, GEMS3KGenerator::f_key_value);

        const std::string legacy_ipm = read_file(data_dir + "/CalcColumn-ipm.dat");
        if (legacy_ipm.find("pa_PSTALL") != std::string::npos)
            throw std::runtime_error("Legacy test IPM input unexpectedly defines pa_PSTALL");

        const std::string marker = "<tMin>";
        const auto marker_position = legacy_ipm.find(marker);
        if (marker_position == std::string::npos)
            throw std::runtime_error("Test IPM input does not contain <tMin>");

        std::string pstall_disabled_ipm = legacy_ipm;
        pstall_disabled_ipm.insert(marker_position, "<pa_PSTALL>  0\n");

        TMultiBase multi;
        std::stringstream first_input(pstall_disabled_ipm);
        multi.read_ipm_format_stream(
            first_input, GEMS3KGenerator::f_key_value, &dch, "pstall_default_test");
        if (multi.base_param()->PSTALL != 0)
        {
            std::cerr << "Explicit pa_PSTALL=0 was not loaded\n";
            return 1;
        }
        multi.multi_kill();

        std::stringstream second_input(legacy_ipm);
        multi.read_ipm_format_stream(
            second_input, GEMS3KGenerator::f_key_value, &dch, "pstall_default_test");
        if (multi.base_param()->PSTALL != 1)
        {
            std::cerr << "Omitted pa_PSTALL retained the previous value\n";
            return 1;
        }

        multi.multi_kill();
        dbr_dch_api::datach_free(&dch);
    }
    catch (const std::exception& error)
    {
        dbr_dch_api::datach_free(&dch);
        std::cerr << error.what() << '\n';
        return 1;
    }

    return 0;
}
