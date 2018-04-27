#include "json.hpp"

#include "rapidjson/document.h"
#include "rapidjson/error/en.h"
#include "rapidjson/writer.h"
#include "rapidjson/prettywriter.h"
#include "rapidjson/stringbuffer.h"
#include <iostream>
#include <sstream>
#include <cstdio>
#include <stdexcept>

JSON_API void read_json(std::string fname, rapidjson::Document &d)
{
    FILE* fp = fopen(fname.c_str(), "rb"); // non-Windows use "r"
    if(!fp)
    {
        std::stringstream str; str << "ERROR [read_json]: file " << fname << " not found!";
        throw std::runtime_error(str.str());
    }

    fseek(fp, 0, SEEK_END);
    size_t length = static_cast<size_t>(ftell(fp));
    //std::cout << "file size = " << length << std::endl;
    fseek(fp, 0, SEEK_SET);
    char* readBuffer = static_cast<char*>(malloc(length + 1));
    size_t readLength = fread(readBuffer, 1, length, fp);
    readBuffer[readLength] = '\0';
    fclose(fp);

    d.Parse<rapidjson::kParseCommentsFlag>(readBuffer);

    if(d.HasParseError())
    {
        std::stringstream str;
        str << "ERROR [read_json]: json document cannot be parsed!" << std::endl;
        str << "\t Parse Error " << d.GetParseError() << " (" << GetParseError_En(d.GetParseError()) << ")" << std::endl;
        str << "\t Error offset = " << d.GetErrorOffset() << '\n';
	free(readBuffer);
        throw std::runtime_error(str.str());
    }    

    if(!d.IsObject())
    {
        std::stringstream str;
        str << "ERROR [read_json]: json document is not an object!";
	free(readBuffer);
        throw std::runtime_error(str.str());
    }
    free(readBuffer);
}

JSON_API bool read_bool(rapidjson::Document const &d, char const *name, bool def)
{
    if(d.HasMember(name))
    {
        if(!d[name].IsBool())
        {
            std::stringstream str;
            str << "ERROR [read_bool]: \"" << name << "\" should be a bool";
            throw std::runtime_error(str.str()); 
        }
        return d[name].GetBool();
    }
    return def;
}

JSON_API int read_int(rapidjson::Document const &d, char const *name, int def)
{
    if(d.HasMember(name))
    {
        if(!d[name].IsInt())
        {
            std::stringstream str;
            str << "ERROR [read_int]: \"" << name << "\" should be an integer";
            throw std::runtime_error(str.str()); 
        }
        return d[name].GetInt();
    }
    return def;
}

JSON_API double read_double(rapidjson::Document const &d, char const *name, double def)
{
    if(d.HasMember(name))
    {
        if(!d[name].IsNumber())
        {
            std::stringstream str;
            str << "ERROR [read_double]: \"" << name << "\" should be a number";
            throw std::runtime_error(str.str()); 
        }
        return d[name].GetDouble();
    }
    return def;
}

JSON_API std::string read_string(rapidjson::Document const &d, char const *name, std::string const &def)
{
    if(d.HasMember(name))
    {
        if(!d[name].IsString())
        {
            std::stringstream str;
            str << "ERROR [read_string]: \"" << name << "\" should be a string";
            throw std::runtime_error(str.str()); 
        }
        return d[name].GetString();
    }
    return def;
}

JSON_API std::vector<double> read_vector_double(
    rapidjson::Document const &d,
    char const *name, 
    std::vector<double> const &def)
{
    if(d.HasMember(name))
    {
        if(!d[name].IsArray())
        {
            std::stringstream str;
            str << "ERROR [read_Vec3d]: \"" << name << "\" should be an array";
            throw std::runtime_error(str.str()); 
        }
        for(rapidjson::SizeType i = 0; i < d[name].Size(); i++)
        {
            if(!d[name][i].IsNumber())
            {
                std::stringstream str;
                str << "ERROR [read_Vec3d]: array \"" << name <<"\" does not contain numbers!";
                throw std::runtime_error(str.str());
            }
        }
        std::vector<double> returned_vector(d[name].Size(),0);
        for(rapidjson::SizeType i = 0; i < d[name].Size(); i++)
            returned_vector[i] = d[name][i].GetDouble();
        return returned_vector;
    }
    return def;
}

JSON_API std::vector<int> read_vector_int(
    rapidjson::Document const &d, 
    char const *name, 
    std::vector<int> const &def)
{
    if(d.HasMember(name))
    {
        if(!d[name].IsArray())
        {
            std::stringstream str;
            str << "ERROR [read_Vec3i]: \"" << name << "\" should be an array";
            throw std::runtime_error(str.str());
        }
        if(d[name].Size()!=3)
        {
            std::stringstream str;
            str << "ERROR [read_Vec3i]: wrong array size for \"" << name <<"\"!";
            throw std::runtime_error(str.str());
        }
        for(rapidjson::SizeType i = 0; i < d[name].Size(); i++)
        {
            if(!d[name][i].IsInt())
            {
                std::stringstream str;
                str << "ERROR [read_Vec3i]: array \"" << name <<"\" does not contain integers!";
                throw std::runtime_error(str.str());
            }
        }
        std::vector<int> returned_vector(d[name].Size(),0);
        for(rapidjson::SizeType i = 0; i < d[name].Size(); i++)
            returned_vector[i] = d[name][i].GetInt();
        return returned_vector;
    }
    return def;
}

JSON_API std::vector<std::string> read_vector_string(
    rapidjson::Document const &d,
    char const *name,
    std::vector<std::string> const &def
)
{
    if(d.HasMember(name))
    {
        if(!d[name].IsArray())
        {
            std::stringstream str;
            str << "ERROR [read_Vec3i]: \"" << name << "\" should be an array";
            throw std::runtime_error(str.str());
        }
        for(rapidjson::SizeType i = 0; i < d[name].Size(); i++)
        {
            if(!d[name][i].IsString())
            {
                std::stringstream str;
                str << "ERROR [read_Vec3i]: array \"" << name <<"\" does not contain strings!";
                throw std::runtime_error(str.str());
            }
        }
        std::vector<std::string> returned_vector(d[name].Size());
        for(rapidjson::SizeType i = 0; i < d[name].Size(); i++)
            returned_vector[i] = d[name][i].GetString();
        return returned_vector;
    }
    return def;
}
