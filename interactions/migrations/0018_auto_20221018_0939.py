# Generated by Django 3.2.15 on 2022-10-18 16:39

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('interactions', '0017_pathway_evidences'),
    ]

    operations = [
        migrations.AddField(
            model_name='ligand',
            name='evidences',
            field=models.TextField(default=''),
        ),
        migrations.AddField(
            model_name='receptor',
            name='evidences',
            field=models.TextField(default=''),
        ),
        migrations.AlterField(
            model_name='pathway',
            name='evidences',
            field=models.TextField(default=''),
        ),
    ]
