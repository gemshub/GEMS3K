#pragma once

#include <vector>
#include <map>
#include <set>
#include <optional>
#include <algorithm>
#include "v_detail.h"

#ifdef USE_NLOHMANNJSON

#include <nlohmann/json.hpp>

class JsonConfigSection
{
public:
    JsonConfigSection(const std::string& filename);

    explicit JsonConfigSection(nlohmann::json arg_json):
        obj_json(arg_json)
    {}

    template <typename T>
    T get_as() const {
        return obj_json.get<T>();
    }

    bool contains(const std::string& name_token ) const
    {
        return obj_json.contains(name_token);
    }

    bool is_string() const
    {
        return obj_json.is_string();
    }

    std::optional<JsonConfigSection> section(const std::string &variable_name_path) const;
    std::vector<JsonConfigSection> to_vector();
    std::string to_string() const
    {
        return obj_json.dump(4);
    }

protected:
    nlohmann::json obj_json;  // Store a pointer to the config file

};

#else


class JsonConfigSectionImpl;

class JsonConfigSection
{
public:
    JsonConfigSection(const std::string& filename);

    explicit JsonConfigSection(const JsonConfigSectionImpl& arg_json);

    template <class T>
    T get_as() const {
        if(std::is_integral_v<T>) {
            return static_cast<T>(as_int());
        }
        else if(std::is_floating_point_v<T>) {
            return static_cast<T>(as_double());
        }
    }

    bool contains(const std::string& name_token ) const;
    bool is_string() const;
    std::optional<JsonConfigSection> section(const std::string &variable_name_path) const;
    std::string to_string() const;
    std::vector<JsonConfigSection> to_vector();

protected:

    // Internal structure of file data
    std::shared_ptr<JsonConfigSectionImpl> impl;

    int64_t as_int() const;
    double as_double() const;

};

template <>  std::string JsonConfigSection::get_as<std::string>() const;
template <>  bool JsonConfigSection::get_as<bool>() const ;
template <>  std::map<std::string, std::vector<std::string>> JsonConfigSection::get_as<std::map<std::string, std::vector<std::string>>>() const;
template <>  std::vector<std::string> JsonConfigSection::get_as<std::vector<std::string>>() const;
template <>  std::set<std::string> JsonConfigSection::get_as<std::set<std::string>>() const;

#endif

/// Definition of json based section settings
class TJsonConfig
{
    JsonConfigSection section_object;

public:

    TJsonConfig(const std::string& filename):
        section_object(filename)
    {}

    explicit TJsonConfig(JsonConfigSection arg_json):
        section_object(arg_json)
    {}

    std::optional<TJsonConfig> section(const std::string &variable_name_path) const
    {
        auto opt_result =  section_object.section(variable_name_path);
        if (!opt_result) {
            return std::nullopt;
        }
        return TJsonConfig(*opt_result);
    }

    /// Check if a section object contains a certain jsonpath.
    bool contains(const std::string& variable_name_path) const
    {
        std::optional<TJsonConfig> opt_result = section(variable_name_path);
        if (!opt_result) {
            return false;
        }
        return true;
    }

    bool is_string() const
    {
        return section_object.is_string();
    }

    template <typename T>
    T get_as() const
    {
        return section_object.get_as<T>();
    }

    template <typename T>
    std::optional<T> value(const std::string& variable_name_path) const
    {
        std::optional<TJsonConfig> opt_result =  section(variable_name_path);
        if (!opt_result) {
            return std::nullopt;
        }
        return opt_result->get_as<T>();
    }

    template <typename T>
    T value_or_default(const std::string& variable_name_path, T default_value) const
    {
        std::optional<T> maybe = value<T>(variable_name_path);
        return maybe.value_or(default_value);
    }

    template <typename T>
    T value_must_exist(const std::string& variable_name_path) const
    {
        std::optional<T> maybe = value<T>(variable_name_path);
        if (!maybe) {
            Error( variable_name_path, "The variable is missing.");
        }
        return *maybe;
    }

    std::vector<TJsonConfig> to_vector()
    {
        auto section_vector = section_object.to_vector();
        std::vector<TJsonConfig> config_as_vector;
        std::transform(
                    section_vector.begin(),
                    section_vector.end(),
                    std::back_inserter(config_as_vector),
                    [](JsonConfigSection item) { return TJsonConfig(item); });
        return config_as_vector;
    }

    /// @brief Dump section to JSON string.
    std::string dump() const
    {
        return  section_object.to_string();
    }

};

/// \class GemsSettings - storing internal preferences
class GemsSettings: public TJsonConfig
{

public:

    /// Task settings file name
    static std::string settings_file_name;

    static std::string logger_section_name;
    static std::set<std::string> default_gems3k_loggers;
    static std::string gems3k_logger_pattern;

    /// Constructor
    explicit GemsSettings(const std::string& config_file_path=GemsSettings::settings_file_name);
    virtual ~GemsSettings();

    /// Remove logging to stdout, logging data only to text file logfile_name
    void gems3k_clear_loggers(const std::string &logfile_name);
    /// Update chemicalfun logger settings
    /// @param use_cout:      show/hide logging to stdout
    ///        logfile_name:  add logging to rotating file name (hide if empty)
    ///        log_level:     set login level for all loggers
    void gems3k_update_loggers(bool use_stdout, const std::string &logfile_name, size_t log_level);

private:

    // logger internal data --------------------------------
    std::optional<TJsonConfig> log_section;
    std::set<std::string> gems3k_loggers = default_gems3k_loggers;
    std::string gems3k_logger_level = "info";

    /// Get the output level for module logger.
    /// @param module_name Name of the logger.
    spdlog::level::level_enum get_level(const std::string& module);

    /// Update/reread spdlog settings
    bool update_logger();
};

GemsSettings& gemsSettings();
