# QueryBuilder 2.5.2 (no jQuery)

This is a modified version of [jQuery QueryBuilder 2.5.2](https://querybuilder.js.org/)
converted to vanilla JavaScript with no external dependencies.

## Original Source

- **Library:** jQuery QueryBuilder 2.5.2
- **Author:** Damien "Mistic" Sorel (http://www.strangeplanet.fr)
- **Original Repository:** https://github.com/mistic100/jQuery-QueryBuilder
- **License:** MIT (https://opensource.org/licenses/MIT)

The original library depended on jQuery, doT.js (template engine), and
jquery-extendext (deep merge utility).

## Modifications

This conversion was performed using AI to remove the jQuery
dependency and strip unused features, while preserving the functionality used by
InterSpec.

### Dependencies removed
- **jQuery** - All DOM manipulation converted to vanilla JavaScript APIs.
- **doT.js** - The 5 HTML templates were converted to plain JavaScript functions.
- **jquery-extendext** - Replaced by a small inline `deepMergeReplace()` utility.

### Features removed (not used by InterSpec)
- All 11 plugins (bt-checkbox, bt-selectpicker, bt-tooltip-errors,
  chosen-selectpicker, filter-description, invert, not-group, sortable,
  sql-support, unique-filter, change-filters)
- `checkbox` and `textarea` input types
- `boolean` and `datetime` type-specific validation
- MomentJS datetime validation
- `optgroup` support for filters and operators
- `multiple` select support
- `value_separator`, `filter.plugin`, `filter.valueGetter`/`valueSetter`
- UMD module wrapper (AMD/CommonJS)

### Features preserved
- Core rule/group creation, deletion, movement
- Filter and operator dropdowns (`text`, `number`, `select`, `radio` inputs)
- Validation with `callback` and `format` (regex)
- `getRules()` / `setRules()` / `validate()` / `reset()` / `clear()` / `destroy()`
- Custom event system (`trigger`/`change`/`on`/`off`)
- Error display (`.has-error` CSS class)
- `string`, `integer`, `double` type coercion
- Condition toggling (AND/OR), group nesting
- Rule/group flags (readonly, no_delete, etc.)

### API changes
- `$('#el').queryBuilder({...})` becomes `QueryBuilder.create(element, {...})`
- `$('#el').queryBuilder('getRules')` becomes `instance.getRules()`
- `$.fn.queryBuilder.constructor.DEFAULTS` becomes `QueryBuilder.DEFAULTS`
- Event binding via `events` option in `QueryBuilder.create()` options

## File Structure

```
js/
  query-builder.js       -- Self-contained vanilla JS, no external dependencies
css/
  query-builder.default.css
  query-builder.default.min.css
  query-builder.dark.css
  query-builder.dark.min.css
  QueryBuilderFakeBootstrap.css
```

CSS files are unchanged from the original QueryBuilder 2.5.2.
