#include <iostream>
#include <fstream>
#include "jsonconfig.h"
#include "v_service.h"
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/rotating_file_sink.h>

#ifdef USE_NLOHMANNJSON

JsonConfigSection::JsonConfigSection(const std::string &filename)
{
    std::string input_str;
    std::fstream f_json(filename, std::ios::in);
    if(f_json.good()) {
    std::stringstream buffer;
    buffer << f_json.rdbuf();
    input_str = buffer.str();
    gems_logger->info( "The configuration file {} has been read.", filename);
    }
    else {
       input_str = "{}";
       gems_logger->debug( "Configuration file {} does not exist.", filename);
    }
    obj_json = nlohmann::json::parse(input_str);
}

std::vector<JsonConfigSection> JsonConfigSection::to_vector()
{
    if (!obj_json.is_array())
    {
        Error( "to_vector", "Must be array type variable.");
    }
    std::vector<JsonConfigSection> json_as_vector;
    std::transform(
                obj_json.begin(),
                obj_json.end(),
                std::back_inserter(json_as_vector),
                [](nlohmann::json item) { return JsonConfigSection(item); });
    return json_as_vector;
}

std::optional<JsonConfigSection> JsonConfigSection::section(const std::string &variable_name_path) const
{
    std::vector<std::string> name_tokens = split(variable_name_path, ".");
    nlohmann::json json_element = obj_json;
    for (const std::string& name_token : name_tokens)
    {
        if (!json_element.contains(name_token))
        {
            return std::nullopt;
        }
        json_element = json_element.at(name_token);
    }
    return JsonConfigSection(json_element);
}

#else


#include <cmath>
#include <string_view>
#include "simdjson/simdjson.h"


std::vector<std::string> as_string_array(const simdjson::dom::element& element);

/// Read fields of structure
class JsonConfigSectionImpl
{
    friend class JsonConfigSection;

public:

    /// Constructor
    explicit JsonConfigSectionImpl( const std::string& json_string )
    {
        try {
            auto input_string = json_string;
            trim(input_string);
            parser = std::make_shared<simdjson::dom::parser>();
            json_element = parser->parse( input_string );
        }
        catch( simdjson::simdjson_error& err )
        {
            gems_logger->error("SimdJson read error : {}.{}", std::to_string(err.error()), err.what());
            Error( std::string("SimdJson read error :") + std::to_string(err.error()), err.what() );
        }

    }

    explicit JsonConfigSectionImpl( std::shared_ptr<simdjson::dom::parser> base_parser,
                                    const simdjson::dom::element& element ):
        parser(base_parser), json_element(element)
    { }

protected:

    std::shared_ptr<simdjson::dom::parser> parser;
    simdjson::dom::element json_element;

};


JsonConfigSection::JsonConfigSection(const std::string &filename)
{
    std::string input_str;
    std::fstream f_json(filename, std::ios::in);
    if(f_json.good()) {
    std::stringstream buffer;
    buffer << f_json.rdbuf();
    input_str = buffer.str();
    gems_logger->info( "The configuration file {} has been read.", filename);
    }
    else {
       input_str = "{}";
       gems_logger->debug("Configuration file {} does not exist.", filename);
    }
    impl = std::make_shared<JsonConfigSectionImpl>(input_str);
}

JsonConfigSection::JsonConfigSection(const JsonConfigSectionImpl &arg_json)
{
    impl = std::make_shared<JsonConfigSectionImpl>(arg_json);
}

bool JsonConfigSection::contains(const std::string &name_token) const
{
    return !impl->json_element.at_key(name_token).error();
}

bool JsonConfigSection::is_string() const
{
    return impl->json_element.is_string();
}

std::string JsonConfigSection::to_string() const
{
    return simdjson::to_string(impl->json_element);
}

int64_t JsonConfigSection::as_int() const
{
    return impl->json_element.get_int64();
}

double JsonConfigSection::as_double() const
{
    return impl->json_element.get_double();
}

std::vector<JsonConfigSection> JsonConfigSection::to_vector()
{
    if (!impl->json_element.is_array()) {
        Error( "to_vector", "Must be array type variable.");
    }
    std::vector<JsonConfigSection> json_as_vector;
    for (const auto&  element : impl->json_element) {
        json_as_vector.push_back(JsonConfigSection(JsonConfigSectionImpl(impl->parser,element)));
    }
    /*std::transform(
                impl->json_element.begin(),
                impl->json_element.end(),
                std::back_inserter(json_as_vector),
                [&](const simdjson::dom::element& item)
    { return JsonConfigSection(JsonConfigSectionImpl(impl->parser,item)); });*/
    return json_as_vector;
}

std::optional<JsonConfigSection> JsonConfigSection::section(const std::string &variable_name_path) const
{
    std::string json_pointer = variable_name_path;
    replaceall( json_pointer, ".", "/");
    json_pointer = "/" + json_pointer;

    auto element = impl->json_element.at_pointer(json_pointer);

    if( element.error() )    {
        return std::nullopt;
    }
    simdjson::dom::element el =impl->json_element.at_pointer(json_pointer);
    return JsonConfigSection(JsonConfigSectionImpl(impl->parser,el));
}

template <>  std::string JsonConfigSection::get_as<std::string>() const {
    std::string_view buf_str = impl->json_element.get_string();
    return std::string( buf_str.begin(), buf_str.end());
}

template <>  bool JsonConfigSection::get_as<bool>() const {
    return impl->json_element.get_bool();
}

std::vector<std::string> as_string_array(const simdjson::dom::element& element)
{
    if (!element.is_array()) {
        Error( "get_as", "Must be array type variable.");
    }
    simdjson::dom::array arr = element.get_array();
    size_t count = arr.size();
    std::vector<std::string> list_data(count);
    size_t index = 0;
    for(std::string_view x : arr) {
        list_data[index++] = x;
    }
    return list_data;
}

template<>
std::map<std::string, std::vector<std::string>> JsonConfigSection::get_as<std::map<std::string, std::vector<std::string>>>() const
{
    std::map<std::string, std::vector<std::string>> data_map;
    if (!impl->json_element.is_object()) {
        Error( "get_as", "Must be object type variable.");
    }
    for (auto field : impl->json_element.get_object()) {
        data_map[std::string(field.key)] = as_string_array(field.value);
    }
    return data_map;
}

template<>
std::vector<std::string> JsonConfigSection::get_as<std::vector<std::string>>() const
{
    return as_string_array(impl->json_element);
}

template<>
std::set<std::string> JsonConfigSection::get_as<std::set<std::string>>() const
{
    auto array = as_string_array(impl->json_element);
    return std::set<std::string>(array.begin(), array.end());
}


#endif

std::string GemsSettings::settings_file_name = "gems3k-config.json";
std::string GemsSettings::logger_section_name = "log";

std::set<std::string> GemsSettings::default_gems3k_loggers = {
    "gems3k", "ipm", "tnode", "kinmet", "solmod",
    #ifdef USE_THERMOFUN
    "chemicalfun", "thermofun"
    #endif
};

std::string GemsSettings::gems3k_logger_pattern("[%n] [%^%l%$] %v");

GemsSettings& gemsSettings()
{
    static  GemsSettings data(GemsSettings::settings_file_name);
    return  data;
}


GemsSettings::GemsSettings(const std::string &config_file_path):
    TJsonConfig(config_file_path)
{
    update_logger();
    //std::fstream f_out( settings_file_name, std::ios::out );
    //f_out << dump();
}

GemsSettings::~GemsSettings()
{
}

spdlog::level::level_enum GemsSettings::get_level(const std::string &module_name)
{
    if (!log_section) {
        return spdlog::level::info;
    }
    auto level_name = log_section->value_or_default<std::string>("module_level."+module_name, gems3k_logger_level);
    return spdlog::level::from_str(level_name);
}

bool GemsSettings::update_logger()
{
    log_section = section(logger_section_name);
    if (!log_section) {
        return false;
    }

    gems3k_loggers = log_section->value_or_default<std::set<std::string>>("modules", {});
    gems3k_loggers.insert(default_gems3k_loggers.begin(), default_gems3k_loggers.end());
    gems3k_logger_level = log_section->value_or_default<std::string>("level", "info");

    std::shared_ptr<spdlog::sinks::stdout_color_sink_mt> stdout_sink;
    std::set<std::string> stdout_module_names;
    auto stdout_section = log_section->section("stdout");
    if(stdout_section)
    {
        stdout_module_names = stdout_section->value_or_default<std::set<std::string>>("modules", {});
        stdout_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
        auto stdout_pattern = stdout_section->value_or_default<std::string>("pattern", gems3k_logger_pattern);
        stdout_sink->set_pattern(stdout_pattern);
    }
    else {
        stdout_module_names = gems3k_loggers;
        stdout_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
        stdout_sink->set_pattern(gems3k_logger_pattern);
    }

    std::shared_ptr<spdlog::sinks::rotating_file_sink_mt> file_sink;
    std::set<std::string> file_module_names;
    auto file_section = log_section->section("file");
    if(file_section)
    {
        file_module_names = file_section->value_or_default<std::set<std::string>>("modules", {});
        auto logfile_path = file_section->value_or_default<std::string>("path", "gems_log.txt");
        auto logfile_size =file_section->value_or_default<size_t>("size", 1048576);
        auto logfile_count = file_section->value_or_default<size_t>("count", 2);
        file_sink = std::make_shared<spdlog::sinks::rotating_file_sink_mt>(
                    logfile_path, logfile_size, logfile_count);
        auto logfile_pattern = file_section->value_or_default<std::string>("pattern", gems3k_logger_pattern);
        file_sink->set_pattern(logfile_pattern);
    }

    for(const auto& lname: gems3k_loggers) {
        auto logger = spdlog::get(lname);
        if (logger) {
            logger->sinks().clear();
            if(stdout_sink && stdout_module_names.find(lname)!=stdout_module_names.end()) {
                logger->sinks().push_back(stdout_sink);
            }
            if(file_sink && file_module_names.find(lname)!=file_module_names.end()) {
                logger->sinks().push_back(file_sink);
            }
            logger->set_level(get_level(lname));
        }
    }

    return true;
}

void GemsSettings::gems3k_update_loggers(bool use_stdout, const std::string& logfile_name, size_t log_level)
{

    spdlog::level::level_enum log_lev = spdlog::level::info;
    if( log_level<7 ) {
        log_lev = static_cast<spdlog::level::level_enum>(log_level);
    }

    auto stdout_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
    stdout_sink->set_pattern(gems3k_logger_pattern);
    std::shared_ptr<spdlog::sinks::rotating_file_sink_mt> file_sink;
    if(!logfile_name.empty()) {
        file_sink = std::make_shared<spdlog::sinks::rotating_file_sink_mt>(logfile_name, 1048576, 1);
        file_sink->set_pattern(gems3k_logger_pattern);
    }

    for(const auto& logger_name: gems3k_loggers) {
        auto logger = spdlog::get(logger_name);
        if(!logger) {
           std::cout <<  logger_name << " logger not connected" << std::endl;
           continue;
        }
        logger->sinks().clear();
        if(use_stdout) {
            logger->sinks().push_back(stdout_sink);
        }
        if(file_sink) {
            logger->sinks().push_back(file_sink);
        }
        logger->set_level(log_lev);
    }
    auto logger = spdlog::get("ipmlog");
    if(logger){ // changed level for file output
       logger->set_level(log_lev);
    }
}

void GemsSettings::gems3k_clear_loggers(const std::string& logfile_name)
{
    std::shared_ptr<spdlog::sinks::rotating_file_sink_mt> file_sink;
    if(!logfile_name.empty()) {
        file_sink = std::make_shared<spdlog::sinks::rotating_file_sink_mt>(logfile_name, 1048576, 1);
    }
    for(const auto& lname: gems3k_loggers) {
        auto logger = spdlog::get(lname);
        if (logger) {
            logger->sinks().clear();
            if(file_sink) {
                logger->sinks().push_back(file_sink);
            }
        }
    }
}
