---
layout: page
title: Publication
description: A publication list with full access (for persnal useage)
keywords: Pulications
comments: true
menu: 链接
permalink: /publication/
---

### Machine learning interatomic potentials.

<ul>
{% for pubs in site.data.pubs %}
  {% if pubs.src == 'ml' %}
  <li><a href="{{ pubs.url }}" target="_blank">{{ pubs.name}}</a></li>
  {% endif %}
{% endfor %}
</ul>


###  Grain boundaries.

<ul>
{% for pubs in site.data.pubs %}
  {% if pubs.src == 'gb' %}
  <li><a href="{{ pubs.url }}" target="_blank">{{ pubs.name}}</a></li>
  {% endif %}
{% endfor %}
</ul>

### Code.

<ul>
{% for pubs in site.data.pubs %}
  {% if pubs.src == 'matsci' %}
  <li><a href="{{ pubs.url }}" target="_blank">{{ pubs.name}}</a></li>
  {% endif %}
{% endfor %}
</ul>
