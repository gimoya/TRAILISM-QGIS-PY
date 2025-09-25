## {{ title }}

{{ summary }}

### Algorithm id
`{{ alg_id }}`

### Inputs
{% for item in inputs %}- **{{ item.label }}** ({{ item.kind }})
{% endfor %}

### Parameters
{% for p in params %}- **{{ p.label }}** ({{ p.kind }}{% if p.default %}, default {{ p.default }}{% endif %})
{% endfor %}

### Outputs
{% for o in outputs %}- **{{ o.label }}**
{% endfor %}

### Usage
```python
import processing
processing.run('{{ alg_id }}', {
  # ...
})
```


