// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#include <odia/MiniJson.h>

#include <cctype>
#include <charconv>
#include <cstddef>
#include <vector>

namespace ODIA::MiniJson
{

  namespace
  {
    constexpr std::size_t npos = std::string_view::npos;

    std::size_t skipSpace(std::string_view s, std::size_t p)
    {
      while (p < s.size() && (s[p] == ' ' || s[p] == '\t' || s[p] == '\n' || s[p] == '\r'))
      {
        ++p;
      }
      return p;
    }

    /// End of the string literal starting at @p p (which must be the opening
    /// quote), one past the closing quote, or npos.
    ///
    /// Escapes matter here rather than being pedantry: a value such as
    /// "bumped by \"schema_version\": 99" contains what a substring scan reads
    /// as a key.
    std::size_t endOfString(std::string_view s, std::size_t p)
    {
      if (p >= s.size() || s[p] != '"') { return npos; }
      for (++p; p < s.size(); ++p)
      {
        if (s[p] == '\\') { ++p; continue; }
        if (s[p] == '"') { return p + 1; }
      }
      return npos;
    }

    std::size_t endOfValue(std::string_view s, std::size_t p);

    /// End of the object or array starting at @p p, one past its closing
    /// bracket, or npos. Nesting and strings are both tracked, so a brace
    /// inside a string does not end the object.
    std::size_t endOfBracketed(std::string_view s, std::size_t p)
    {
      const char open = s[p];
      const char close = open == '{' ? '}' : ']';
      int depth = 0;
      for (; p < s.size(); ++p)
      {
        if (s[p] == '"')
        {
          const auto e = endOfString(s, p);
          if (e == npos) { return npos; }
          p = e - 1;
          continue;
        }
        if (s[p] == open) { ++depth; }
        else if (s[p] == close)
        {
          --depth;
          if (depth == 0) { return p + 1; }
        }
      }
      return npos;
    }

    std::size_t endOfValue(std::string_view s, std::size_t p)
    {
      p = skipSpace(s, p);
      if (p >= s.size()) { return npos; }
      if (s[p] == '"') { return endOfString(s, p); }
      if (s[p] == '{' || s[p] == '[') { return endOfBracketed(s, p); }
      const auto e = s.find_first_of(",}]\n\t\r ", p);
      return e == npos ? s.size() : e;
    }

    /// Position of the value of member @p key of the object starting at @p p,
    /// or npos. Only that object's own members are considered.
    std::size_t member(std::string_view s, std::size_t p, std::string_view key)
    {
      p = skipSpace(s, p);
      if (p >= s.size() || s[p] != '{') { return npos; }
      ++p;
      while (true)
      {
        p = skipSpace(s, p);
        if (p >= s.size() || s[p] == '}') { return npos; }
        if (s[p] != '"') { return npos; }
        const auto name_end = endOfString(s, p);
        if (name_end == npos) { return npos; }
        const auto name = s.substr(p + 1, name_end - p - 2);

        p = skipSpace(s, name_end);
        if (p >= s.size() || s[p] != ':') { return npos; }
        const auto value = skipSpace(s, p + 1);
        if (name == key) { return value; }

        const auto value_end = endOfValue(s, value);
        if (value_end == npos) { return npos; }
        p = skipSpace(s, value_end);
        if (p < s.size() && s[p] == ',') { ++p; }
      }
    }

    std::string unescape(std::string_view raw)
    {
      std::string out;
      out.reserve(raw.size());
      for (std::size_t i = 0; i < raw.size(); ++i)
      {
        if (raw[i] != '\\' || i + 1 >= raw.size()) { out.push_back(raw[i]); continue; }
        switch (raw[++i])
        {
          case 'n': out.push_back('\n'); break;
          case 't': out.push_back('\t'); break;
          case 'r': out.push_back('\r'); break;
          case 'b': out.push_back('\b'); break;
          case 'f': out.push_back('\f'); break;
          // \uXXXX is passed through as written. Nothing this reads from a
          // census or a version string needs it, and a wrong decoding would be
          // worse than a visible escape.
          default:  out.push_back(raw[i]); break;
        }
      }
      return out;
    }
  } // namespace

  std::optional<std::string> get(std::string_view json,
                                 std::initializer_list<std::string_view> path)
  {
    std::size_t p = skipSpace(json, 0);
    for (const auto& key : path)
    {
      p = member(json, p, key);
      if (p == npos) { return std::nullopt; }
    }
    p = skipSpace(json, p);
    if (p >= json.size()) { return std::nullopt; }
    if (json[p] == '{' || json[p] == '[') { return std::nullopt; }
    if (json[p] == '"')
    {
      const auto e = endOfString(json, p);
      if (e == npos) { return std::nullopt; }
      return unescape(json.substr(p + 1, e - p - 2));
    }
    const auto e = endOfValue(json, p);
    if (e == npos || e <= p) { return std::nullopt; }
    return std::string(json.substr(p, e - p));
  }

  std::optional<long long> getInt(std::string_view json,
                                  std::initializer_list<std::string_view> path)
  {
    const auto raw = get(json, path);
    if (!raw) { return std::nullopt; }
    long long value = 0;
    const auto* begin = raw->data();
    const auto* end = begin + raw->size();
    const auto result = std::from_chars(begin, end, value);
    // The whole token must be the number: "1x" is not 1, and neither is "1.5".
    if (result.ec != std::errc{} || result.ptr != end) { return std::nullopt; }
    return value;
  }

  std::optional<std::string> getString(std::string_view json,
                                       std::initializer_list<std::string_view> path)
  {
    return get(json, path);
  }

} // namespace ODIA::MiniJson
