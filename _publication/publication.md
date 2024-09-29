---
layout: page
title: Publication
description: A publication list with full access (for persnal useage)
keywords: Pulications
comments: true
menu: 链接
permalink: /pubs/
---

> Below are some links I considered to be useful. 

<ul>
{% for pubs in site.data.pubs %}
  {% if pubs.src == 'science' %}
  <li><a href="{{ link.url }}" target="_blank">{{ link.name}}</a></li>
  {% endif %}
{% endfor %}
</ul>