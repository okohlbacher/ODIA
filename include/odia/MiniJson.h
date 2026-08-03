// Copyright (c) 2026, Oliver Kohlbacher and the ODIA authors.
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <initializer_list>
#include <optional>
#include <string>
#include <string_view>

namespace ODIA
{

  /// Just enough JSON to read a scalar at a known path.
  ///
  /// This replaces a substring scan that looked for `"key"` and took the next
  /// token after the colon. That scan was not scoped to anything: asked for
  /// `counts.transitions.total` it found the first `"transitions"` anywhere in
  /// the document and then the first `"total"` after it, at any depth. On a
  /// real `library/metadata.json`, whose per-table blocks carry `target` and
  /// `decoy` but not always `total`, it walked past the block it was asked
  /// about and into a later one -- so a bundle declaring 8 transitions and
  /// carrying 4 reported that its census agreed. That check exists to catch a
  /// truncated read before an hour of extraction; a scan that can pass a
  /// truncated file is worse than no check.
  ///
  /// Deliberately not a JSON library: it validates nothing it is not asked
  /// about and builds no tree. It knows string escapes and nesting, which is
  /// exactly what the substring scan did not.
  namespace MiniJson
  {
    /// The scalar at @p path, as it appears in the document, or nullopt if the
    /// path does not resolve to one. Strings are returned unquoted and
    /// unescaped only for the simple escapes; a value that is an object or an
    /// array is not a scalar and yields nullopt.
    std::optional<std::string> get(std::string_view json,
                                   std::initializer_list<std::string_view> path);

    /// Convenience wrappers; nullopt also when the scalar does not parse.
    std::optional<long long> getInt(std::string_view json,
                                    std::initializer_list<std::string_view> path);
    std::optional<std::string> getString(std::string_view json,
                                         std::initializer_list<std::string_view> path);
  }

} // namespace ODIA
