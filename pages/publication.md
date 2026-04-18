---
layout: page
title: Publication
description: A publication list with full access (for persnal useage)
keywords: Pulications
comments: true
menu: 链接
permalink: /publication/
---

---
layout: page
title: Publication
---

## Machine learning interatomic potentials

{% bibliography -f papers --query @*[keywords=ml] %}


## All publications
<ul>
{% for pubs in site.data.pubs %}
  <li><a href="{{ pubs.url }}" target="_blank">{{ pubs.name }}</a></li>
{% endfor %}
</ul>