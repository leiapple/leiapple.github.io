---
layout: page
title: Links
description: A collection of useful links
keywords: 友情链接
comments: true
menu: 链接
permalink: /links/
---

> Adademic groups in computational materials.

<ul>
{% for link in site.data.links %}
  {% if link.src == 'science' %}
  <li><a href="{{ link.url }}" target="_blank">{{ link.name}}</a></li>
  {% endif %}
{% endfor %}
</ul>


> Courses/Resources in Machine learning interatomic potentials.

<ul>
{% for link in site.data.links %}
  {% if link.src == 'course' %}
  <li><a href="{{ link.url }}" target="_blank">{{ link.name}}</a></li>
  {% endif %}
{% endfor %}
</ul>