# Generated by Django 4.0.2 on 2022-08-18 22:46

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('interactions', '0008_cellclas_cell_type'),
    ]

    operations = [
        migrations.AddField(
            model_name='pathwayandcelltype',
            name='averageScore',
            field=models.FloatField(default=0),
            preserve_default=False,
        ),
    ]
