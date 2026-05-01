# {{ options.spectrum_title }}

**Status:** {{ status.fail_reason }}{% if not status.success %} — {{ status.error_message }}{% endif %}{{ "" }}
**χ²/dof:** {{ chi2_str }} / {{ dof }} = {{ chi2_per_dof_str }}
**Solver:** {{ timing.duration_str }} ({{ timing.solve_calls }} fcn evals + {{ timing.cov_calls }} for covariance)

<h2>Relative-Efficiency equations</h2>
{% for curve in rel_eff_curves %}
- **{{ curve.name }}** ({{ curve.rel_eff_eqn_type }}): `{{ curve.equation_text }}`
{% endfor %}

<h2>Activities</h2>
{% for curve_acts in relative_activities %}
<h3>Curve {{ curve_acts.curve_index }}: {{ length(curve_acts.nuclides) }} nuclide(s)</h3>

| Nuclide | Rel. Activity | Mass % | Enrichment |
|---|---|---|---|
{% for nuc in curve_acts.nuclides %}
| {{ nuc.name }} | {{ printCompact(nuc.rel_activity, 6) }} ± {{ printCompact(nuc.rel_activity_uncertainty, 6) }} | {% if existsIn(nuc, "total_mass_fraction") %}{{ pct(nuc.total_mass_fraction, 4) }}% {% endif %}| {% if existsIn(nuc, "enrichment") %}{{ pct(nuc.enrichment, 4) }}% {% endif %}|
{% endfor %}
{% endfor %}

{% if length(warnings) > 0 %}
<h2>Warnings</h2>
{% for w in warnings %}
- {{ w }}
{% endfor %}
{% endif %}
