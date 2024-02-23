#pragma once
#include "Config.hpp"
#include <yaml-cpp/yaml.h>
#include "Geometry.hpp"

class InputParser
{

    YAML::Node root_node;
    string input_filename;

public:
    InputParser(const string &input_filename);

    void parse_config(Config &config) const;

    void create_geometry(unique_ptr<Geometry> &geometry) const;

private:
    void parse_yaml_config_options(Config &config) const;

    string option_not_specified_msg(string root_name, string option_name) const;

    template <typename T>
    T read_required_option(string root_name, string option_name) const
    {
        if (root_node[root_name][option_name])
        {
            return root_node[root_name][option_name].as<T>();
        }
        else
        {
            throw std::runtime_error(option_not_specified_msg(root_name, option_name));
        }
    }

    template <typename T>
    T read_optional_option(string root_name, string option_name, T default_value) const
    {
        if (root_node[root_name][option_name])
        {
            return root_node[root_name][option_name].as<T>();
        }
        else
        {
            return default_value;
        }
    }

    template <typename EnumType>
    static EnumType lookup_enum_option_map(const std::map<string, EnumType> &map,
                                           const string &key,
                                           string option_name,
                                           const string &option_parent_name = "")
    {
        if (map.count(key) == 1)
            return map.at(key);
        else
        {
            string keys;
            for (const auto &pair : map)
                keys += "'" + pair.first + "'\n";
            if (option_parent_name.size() > 0)
                option_name = option_parent_name + ": " + option_name;
            throw std::runtime_error("Illegal value '" + key + "' specified for setting '" +
                                     option_name + "'. Legal values are:\n" + keys);
        }
    }

    template <typename EnumType>
    EnumType read_required_enum_option(string root_name, string option_name, const std::map<string, EnumType> &enum_map) const
    {
        if (root_node[root_name][option_name])
        {
            return lookup_enum_option_map(enum_map, root_node[root_name][option_name].as<string>(), option_name);
        }
        else
        {
            throw std::runtime_error(option_not_specified_msg(root_name, option_name));
        }
    }

    template <typename EnumType>
    EnumType read_optional_enum_option(string root_name, string option_name, const std::map<string, EnumType> &enum_map, EnumType default_value) const
    {
        if (root_node[root_name][option_name])
        {
            return lookup_enum_option_map(enum_map, root_node[root_name][option_name].as<string>(), option_name);
        }
        else
        {
            return default_value;
        }
    }
};